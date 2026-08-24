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

from aleph.laws import units as U
from aleph.laws.fiber_network import FiberNetwork, build_fiber_network


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
    cortex_thickness_um: float = 0.0        # radial shell thickness [µm]. 0.0 = the historical ZERO-thickness
                                            # shell (all filaments at exactly R → crossing arcs' nodes coincide
                                            # in radius → ~64k WCA interpenetrations, the build artifact that
                                            # dominates the resting-convergence residual). >0 disperses each
                                            # filament onto its own radius in [R−t/2, R+t/2] so crossing
                                            # filaments sit at DIFFERENT radii and no longer interpenetrate.
                                            # Physical: h_cortex ≈ 0.2 µm (KB-3.1/3.5) ≫ σ_EV = 7 nm. The 2-D
                                            # shell was the simplification; the real cortex has thickness.

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


def _resolve_cross_fiber_overlaps(pos: np.ndarray, fiber_of: np.ndarray, sigma_ev_um: float,
                                  max_sweeps: int, relax: float = 1.0) -> tuple[np.ndarray, int, float]:
    """Push CROSS-FILAMENT interpenetrating nodes apart to the WCA contact distance (build-time relaxation).

    A soft-sphere separation: each sweep moves every different-filament pair closer than ``r_c = 2^(1/6)·σ_EV``
    apart by half their overlap along the pair-separation direction (both nodes symmetrically). Same-filament
    neighbours are never pushed (they are the bonded arc). Converges in a handful of sweeps because the cortex
    now has radial thickness to settle into; the residual node motion is sub-σ_EV (≪ the 0.5 µm mesh), so the
    fiber geometry / on-shell property is preserved and every node stays within the physical cortex shell. The
    relaxation target is the WCA cutoff — a param-free geometric criterion, NOT tuned to a force gate.

    Returns ``(relaxed_pos, sweeps_used, max_node_move_um)``.
    """
    from scipy.spatial import cKDTree
    r_c = 2.0 ** (1.0 / 6.0) * float(sigma_ev_um)
    r_target = r_c * 1.01                                       # clear the cutoff (settle just past it, not on it)
    p = np.ascontiguousarray(pos, np.float64).copy()
    sweeps = 0
    for sweeps in range(1, int(max_sweeps) + 1):
        pairs = cKDTree(p).query_pairs(r_c, output_type="ndarray")
        if pairs.shape[0]:
            pairs = pairs[fiber_of[pairs[:, 0]] != fiber_of[pairs[:, 1]]]
        if pairs.shape[0] == 0:
            break
        i, j = pairs[:, 0], pairs[:, 1]
        d = p[i] - p[j]
        r = np.maximum(np.linalg.norm(d, axis=1), 1e-9)
        push = (float(relax) * 0.5 * (r_target - r) / r)[:, None] * d
        disp = np.zeros_like(p)
        np.add.at(disp, i, push)
        np.add.at(disp, j, -push)
        p = p + disp
    return p, sweeps, float(np.linalg.norm(p - pos, axis=1).max())


def _cosine_span_weights(span: int) -> np.ndarray:
    r"""Symmetric smooth taper ``w[0..span]`` for the radial bump: ``w[0]=1`` (peak at the crossing node),
    tapering as ``cos²(½π·k/(span+1))`` so ``w`` reaches ~0 one node past the window edge. A cos² (Hann-like)
    window has zero slope at the peak, so the bump adds NO curvature spike at the crossing — the whole point of
    spreading the separation over neighbours instead of shoving a single node."""
    k = np.arange(int(span) + 1, dtype=np.float64)
    return np.cos(0.5 * np.pi * k / (int(span) + 1)) ** 2


def _resolve_cross_fiber_overlaps_radial_span(
    pos: np.ndarray, fiber_offsets: np.ndarray, centre: np.ndarray, sigma_ev_um: float,
    max_sweeps: int, span: int = 2, relax: float = 1.0,
    anchor_pairs: np.ndarray | None = None) -> tuple[np.ndarray, int, float]:
    r"""SMOOTHNESS-PRESERVING overlap resolution: separate crossing fibers OUT-OF-PLANE (sphere-radial) with the
    displacement spread over a span of neighbouring nodes, so neither fiber is kinked (opt-in alternative to the
    in-plane transverse push in :func:`_resolve_cross_fiber_overlaps`).

    The transverse push moves a single crossing node ~½·overlap IN the shell tangent plane — off its smooth
    great-circle arc — creating a per-node turn-angle kink whose cytosim bending force ``F=α·d`` (α=κ/seg³)
    blows up ×288 at the fine 75 nm mesh (the fine-cortex GATE-A residual plateau). This mode instead:

      1. resolves each sub-``r_c`` cross-fiber pair by displacing the two nodes along their LOCAL RADIAL
         direction (normal to the shell) — the outer-radius node further out, the inner-radius node further in —
         so the fibers separate in 3D while each node's IN-PLANE (arc) position is untouched; and
      2. spreads every crossing node's radial displacement over a ``±span``-node cosine bump ALONG its own fiber
         (peak at the crossing, tapering smoothly to ~0 at the edge), so the arc bends as a gentle low-curvature
         bump rather than a single-node spike. Per-node turn angle ≈ bump-amplitude / (span·seg), which the span
         makes small even at fine seg, instead of the transverse ≈ amplitude / seg.

    Same param-free geometric target as the transverse mode (settle just past ``r_c = 2^(1/6)·σ_EV``), iterated
    to convergence. ``anchor_pairs`` (Arp2/3 branch welds) are excluded from the push and re-welded each sweep,
    mirroring :func:`ac.weave.regions._relax_arp23_overlaps`. Returns ``(relaxed_pos, sweeps_used, max_move_um)``.
    """
    from scipy.spatial import cKDTree
    r_c = 2.0 ** (1.0 / 6.0) * float(sigma_ev_um)
    r_target = r_c * 1.01
    p = np.ascontiguousarray(pos, np.float64).copy()
    c = np.asarray(centre, np.float64).reshape(3)
    off = np.asarray(fiber_offsets, np.int64)
    F = len(off) - 1
    lens = np.diff(off)
    fiber_of = np.repeat(np.arange(F), lens)                       # per-node filament id
    N = p.shape[0]
    gidx = np.arange(N)
    node_in_fiber = (gidx - off[fiber_of]).astype(np.int64)        # index within own fiber [0, len-1]
    fiber_len = lens[fiber_of]                                     # own-fiber node count, per node
    w = _cosine_span_weights(span)
    span = int(span)
    if anchor_pairs is not None and anchor_pairs.size:
        d0 = anchor_pairs[:, 0].astype(np.int64)
        bn = anchor_pairs[:, 1].astype(np.int64)
        anchor_key = {(min(int(a), int(b)), max(int(a), int(b))) for a, b in anchor_pairs}
    else:
        d0 = bn = np.zeros(0, np.int64)
        anchor_key = set()
    sweeps = 0
    for sweeps in range(1, int(max_sweeps) + 1):
        pairs = cKDTree(p).query_pairs(r_c, output_type="ndarray")
        if pairs.shape[0]:
            pairs = pairs[fiber_of[pairs[:, 0]] != fiber_of[pairs[:, 1]]]
        if pairs.shape[0] and anchor_key:
            keep = np.fromiter(((int(a), int(b)) not in anchor_key for a, b in pairs), bool, pairs.shape[0])
            pairs = pairs[keep]
        if pairs.shape[0] == 0:
            break
        i, j = pairs[:, 0], pairs[:, 1]
        d = p[i] - p[j]
        r = np.maximum(np.linalg.norm(d, axis=1), 1e-9)
        need = float(relax) * 0.5 * (r_target - r)                # per-node radial shove (half each side)
        rad = np.linalg.norm(p - c, axis=1)
        sign_i = np.where(rad[i] >= rad[j], 1.0, -1.0)            # outer node → +radial, inner → −radial
        amp = np.zeros(N)
        np.add.at(amp, i, sign_i * need)
        np.add.at(amp, j, -sign_i * need)
        u = (p - c) / np.maximum(np.linalg.norm(p - c, axis=1, keepdims=True), 1e-12)   # radial unit / node
        disp = np.zeros_like(p)
        for k in range(-span, span + 1):                          # spread each node's amp onto its arc-neighbours
            wk = w[abs(k)]
            if wk == 0.0:
                continue
            valid = (node_in_fiber + k >= 0) & (node_in_fiber + k <= fiber_len - 1)
            src = gidx[valid]
            dst = src + k                                         # same fiber ⇒ neighbour is a contiguous index
            disp[dst] += (wk * amp[src])[:, None] * u[dst]        # apply along the NEIGHBOUR's radial direction
        p = p + disp
        if d0.size:                                               # re-weld daughter base ↔ branch vertex
            p[d0] = p[bn]
    return p, sweeps, float(np.linalg.norm(p - pos, axis=1).max())


def overlap_free_shell(pos: np.ndarray, fiber_offsets: np.ndarray, centre: np.ndarray,
                       thickness_um: float, sigma_ev_um: float = 0.007, *,
                       rng: np.random.Generator | None = None, max_sweeps: int = 60,
                       overlap_mode: str = "transverse", overlap_span: int = 2) -> tuple[np.ndarray, dict]:
    """Make ANY already-placed fiber shell overlap-free: radially disperse across ``thickness_um`` then relax.

    Post-placement counterpart of the in-builder path — usable on a shell built elsewhere (e.g. the ac/ woven
    cortex). Each fiber is first scaled to its own radius in ``[R−t, R]`` about ``centre`` (INWARD only: keeps
    the outer surface at R so membrane↔cortex ERM tethers stay radial/non-degenerate; breaks the crossing-node
    radial coincidence), then :func:`_resolve_cross_fiber_overlaps` clears the residual. Returns
    ``(new_pos, stats)``. ``thickness_um=0`` disperses nothing and only relaxes (rarely enough on a 2-D shell).
    """
    if rng is None:
        rng = np.random.default_rng(0)
    p = np.ascontiguousarray(pos, np.float64).copy()
    c = np.asarray(centre, np.float64)
    off = np.asarray(fiber_offsets, np.int64)
    F = len(off) - 1
    fiber_of = np.repeat(np.arange(F), np.diff(off))
    if thickness_um > 0.0:
        r = np.linalg.norm(p - c, axis=1, keepdims=True)
        scale = np.ones(p.shape[0])
        dR = -rng.random(F) * float(thickness_um)                        # per-fiber INWARD radial offset [−t, 0]
        rad = np.maximum(r[:, 0], 1e-9)
        scale = (rad + dR[fiber_of]) / rad                              # scale each node to its fiber's radius
        p = c + (p - c) * scale[:, None]
    if overlap_mode == "transverse":
        p, sweeps, move = _resolve_cross_fiber_overlaps(p, fiber_of, sigma_ev_um, max_sweeps)
    elif overlap_mode == "radial_span":
        p, sweeps, move = _resolve_cross_fiber_overlaps_radial_span(
            p, off, c, sigma_ev_um, max_sweeps, span=overlap_span)
    else:
        raise ValueError(f"overlap_mode must be 'transverse' or 'radial_span', got {overlap_mode!r}")
    return p, {"overlap_relax_sweeps": int(sweeps), "overlap_relax_max_move_um": float(move),
               "cortex_thickness_um": float(thickness_um), "overlap_mode": overlap_mode}


def count_interpenetrations(net: FiberNetwork, sigma_ev_um: float = 0.007) -> dict:
    """Count CROSS-FILAMENT node pairs closer than the WCA contact distance 2^(1/6)·σ_EV (build diagnostic).

    Same-filament neighbours (bonded/adjacent) are excluded — only DIFFERENT filaments interpenetrating counts,
    because that is the build artifact the steric force fights at t0. Uses a KD-tree (cheap vs the O(N²) naive
    count); returns the number of interpenetrating nodes and the offending pair count. σ_EV = 7 nm (F-actin).
    """
    from scipy.spatial import cKDTree
    r_contact = 2.0 ** (1.0 / 6.0) * float(sigma_ev_um)
    pos = np.ascontiguousarray(net.pos, np.float64)
    off = np.asarray(net.fiber_offsets, np.int64)
    fiber_of = np.repeat(np.arange(len(off) - 1), np.diff(off))     # per-node filament id
    pairs = cKDTree(pos).query_pairs(r_contact, output_type="ndarray")
    if pairs.shape[0] == 0:
        return {"n_interpenetrating_nodes": 0, "n_cross_pairs": 0, "r_contact_um": r_contact}
    cross = fiber_of[pairs[:, 0]] != fiber_of[pairs[:, 1]]          # keep only different-filament pairs
    cross_pairs = pairs[cross]
    interp_nodes = np.unique(cross_pairs.ravel()).size if cross_pairs.shape[0] else 0
    return {"n_interpenetrating_nodes": int(interp_nodes), "n_cross_pairs": int(cross_pairs.shape[0]),
            "r_contact_um": r_contact}


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
                         nematic_S: float = 1.0,
                         length_dist: str = "mono",
                         resolve_overlaps: bool = False,
                         sigma_ev_um: float = 0.007,
                         max_relax_sweeps: int = 60,
                         overlap_mode: str = "transverse",
                         overlap_span: int = 2) -> tuple[FiberNetwork, dict]:
    """Assemble the cortex fiber network on the sphere (FF units, µm).

    Args:
        params: cortex parameters (defaults = H.3 production cortex).
        rng: random generator (default seed 0 for reproducibility).
        n_filaments: optional override of ``params.n_filaments`` (e.g. a small prototype).
        length_dist: per-filament contour-length model.
            ``"mono"`` (default) — every filament is exactly ``L_filament_um`` (the historical,
            γ-validated cortex; kept as the default so validation is untouched).
            ``"exponential"`` — draw each filament length from an exponential of the SAME mean,
            clipped to the KB-3.18 cortical range 1–10 µm. Rationale: real cortical F-actin is
            length-DISTRIBUTED, not monodisperse (KB-3.18: linear F-actin 1–10 µm; Fritzsche 2013:
            two sub-populations of different mean length), and the steady-state length distribution
            of filaments under stochastic capping/severing IS exponential (Edelstein-Keshet;
            Mogilner-Oster). Mean-preserving → the γ calibrated at the mean length is unchanged in
            expectation; only the physical spread is added.

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
    if length_dist == "exponential":
        # KB-3.18: cortical F-actin length is DISTRIBUTED (1–10 µm), not a single value.
        L_f = np.clip(rng.exponential(params.L_filament_um, F), seg, 10.0)  # mean-preserving
        nb_f = np.maximum(2, np.rint(L_f / seg).astype(int) + 1)            # beads = L/ℓ₀ + 1, ≥2
    elif length_dist == "mono":
        nb_f = np.full(F, nb, dtype=int)
    else:
        raise ValueError(f"length_dist must be 'mono' or 'exponential', got {length_dist!r}")
    # Per-filament radius: disperse across the physical cortex thickness so crossing filaments do not share a
    # radius (removes the WCA build-interpenetration artifact). t=0 reproduces the historical zero-thickness
    # shell exactly (R_f ≡ R). The distribution is uniform INWARD across [R−t, R]: the cortex OUTER surface
    # stays at R (= R_CORTEX_UM, one submembranous gap inside the membrane) and the physical thickness extends
    # toward the cell interior. A symmetric ±t/2 spread would push the outer nodes out to R+t/2 = the membrane
    # radius, collapsing the membrane↔cortex ERM tethers to ~zero length / tangential direction so they cannot
    # transmit the radial turgor load (2026-07-22 resting-baseline geometry diagnosis). Mesh-scale, param-free.
    t_shell = float(params.cortex_thickness_um)
    R_f = np.full(F, R) if t_shell <= 0.0 else R - rng.random(F) * t_shell
    fibers = [_great_circle_arc(com_dirs[f], tangents[f], float(R_f[f]), int(nb_f[f]), seg) for f in range(F)]

    # Overlap-free build (opt-in): disperse across the shell thickness already dropped crossing-node coincidence
    # ~18×; a short soft-sphere relaxation clears the residual so the cortex starts with ZERO WCA steric force.
    relax_sweeps, relax_move_um = 0, 0.0
    if resolve_overlaps:
        counts = np.array([a.shape[0] for a in fibers], np.int64)
        fiber_of = np.repeat(np.arange(F), counts)
        pos_flat = np.concatenate(fibers, axis=0)
        edges = np.concatenate([[0], np.cumsum(counts)])
        if overlap_mode == "transverse":
            pos_flat, relax_sweeps, relax_move_um = _resolve_cross_fiber_overlaps(
                pos_flat, fiber_of, sigma_ev_um, max_relax_sweeps)
        elif overlap_mode == "radial_span":
            # smoothness-preserving out-of-plane resolution (cortex sphere is centred at the origin)
            pos_flat, relax_sweeps, relax_move_um = _resolve_cross_fiber_overlaps_radial_span(
                pos_flat, edges.astype(np.int64), np.zeros(3), sigma_ev_um, max_relax_sweeps, span=overlap_span)
        else:
            raise ValueError(f"overlap_mode must be 'transverse' or 'radial_span', got {overlap_mode!r}")
        fibers = [pos_flat[edges[f]:edges[f + 1]] for f in range(F)]

    net = build_fiber_network(fibers, kappa=params.kappa)

    radii = np.linalg.norm(net.pos, axis=1)
    beads_f = np.diff(net.fiber_offsets)                       # actual beads per fiber
    Lf = (beads_f - 1) * seg                                    # per-fiber contour length [µm]
    meta = {
        "n_filaments": F,
        "n_nodes": net.n_nodes,
        "length_dist": length_dist,
        "L_fil_mean_um": float(Lf.mean()), "L_fil_std_um": float(Lf.std()),
        "L_fil_min_um": float(Lf.min()), "L_fil_max_um": float(Lf.max()),
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
        "cortex_thickness_um": t_shell,
        "resolve_overlaps": resolve_overlaps,
        "overlap_relax_sweeps": int(relax_sweeps),
        "overlap_relax_max_move_um": float(relax_move_um),
    }
    if resolve_overlaps:
        meta["n_interpenetrating_nodes"] = count_interpenetrations(net, sigma_ev_um)["n_interpenetrating_nodes"]
    return net, meta


def equatorial_circumference_um(params: CortexParams) -> float:
    """Equatorial circumference 2πR [µm] — the denominator for the method-of-planes γ (Stage 6d)."""
    return 2.0 * np.pi * params.R_um
