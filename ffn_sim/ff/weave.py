"""Unified actin-architecture builder (task a, increment 1) — ONE ``weave()`` for every structure.

``weave(spec, rng)`` takes an :class:`~ff.architecture_spec.ArchitectureSpec` and builds the woven
network: (1) place filaments on the manifold per the filament spec, (2) pair crosslinkers (KDTree
near-pairs on DIFFERENT filaments, filtered by the crosslinker bind-mode angle), (3) place motors.
The SAME function produces the isotropic cortical shell (CORTEX, manifold="sphere") and a parallel
bundle (FILOPODIUM, manifold="bundle") — no structure-specific branch beyond the placement + the
bind-mode angle filter. Returns a :class:`~ff.gamma_floor.CrosslinkedCortex`, so the existing FF relax
/ γ / metrics consume the result unchanged.

The "sphere" path mirrors ``gamma_floor.build_crosslinked_cortex`` (same mesoscale reach + α-actinin/
filamin split + KDTree pairing) so ``weave(CORTEX)`` reproduces the H.3 γ-floor cortex (statistical
parity, `tests/ff/test_weave.py`). Increment 1 scope: CORTEX + FILOPODIUM only (PI-gated table).
"""

from __future__ import annotations

import numpy as np

from ffn_sim.ff.architecture_spec import ArchitectureSpec
from ffn_sim.ff.cortex_assembly import CortexParams, _random_unit_vectors, build_cortex_network
from ffn_sim.ff.fiber_network import build_fiber_network
from ffn_sim.ff.gamma_floor import (
    ALPHA_ACTININ,
    FILAMIN,
    CrosslinkedCortex,
    _cross_fiber_pairs,
    _node_fiber_map,
    mesoscale_reach,
)


def _fiber_tangents(net, pairs: np.ndarray) -> np.ndarray:
    """Unit tangent of each node's fiber (from its segment), for the bind-mode angle filter."""
    fib = _node_fiber_map(net)
    off = net.fiber_offsets
    t = np.zeros((net.n_nodes, 3))
    for f in range(net.n_fibers):
        a, b = int(off[f]), int(off[f + 1])
        seg = net.pos[a + 1:b] - net.pos[a:b - 1]
        # node tangent = mean of adjacent segment directions
        d = seg / (np.linalg.norm(seg, axis=1, keepdims=True) + 1e-12)
        tn = np.zeros((b - a, 3)); tn[:-1] += d; tn[1:] += d
        t[a:b] = tn / (np.linalg.norm(tn, axis=1, keepdims=True) + 1e-12)
    return t


def _filter_bind_mode(net, pairs: np.ndarray, bind_mode: str, theta_max_deg: float = 30.0) -> np.ndarray:
    """Keep crosslink pairs whose two filament tangents satisfy the bind mode (parallel / perp / any)."""
    if bind_mode == "any" or pairs.shape[0] == 0:
        return pairs
    t = _fiber_tangents(net, pairs)
    cos = np.abs(np.einsum("ij,ij->i", t[pairs[:, 0]], t[pairs[:, 1]]))
    cmax = np.cos(np.deg2rad(theta_max_deg))
    keep = cos >= cmax if bind_mode == "parallel" else cos <= np.cos(np.deg2rad(90 - theta_max_deg))
    return pairs[keep]


def _build_bundle(spec: ArchitectureSpec, rng: np.random.Generator):
    """Parallel-bundle placement (filopodium/microvillus core): N straight filaments along +z, packed
    on a 2D hexagonal-ish cross-section at the lit crosslinker spacing, uniform polarity."""
    fs = spec.filament
    nb = int(round(fs.length_um / fs.seg_um)) + 1
    spacing_um = fs.bundle_spacing_um if fs.bundle_spacing_um is not None else 0.008  # ~8 nm fascin / ~12 nm SF-MV
    n = fs.n_filaments
    # hexagonal-ish packing of n points in a 2D disk cross-section
    k = int(np.ceil(np.sqrt(n)))
    xs, ys = [], []
    for i in range(k):
        for j in range(k):
            xs.append((i + 0.5 * (j % 2)) * spacing_um); ys.append(j * spacing_um * np.sqrt(3) / 2)
    xs, ys = np.array(xs[:n]), np.array(ys[:n])
    xs -= xs.mean(); ys -= ys.mean()
    z = (np.arange(nb) - (nb - 1) / 2.0) * fs.seg_um
    # GRADED antiparallel polarity for "mixed" (stress fiber): half the filaments run +z, half −z, so
    # barbed/+ ends cluster at the two ends and pointed/− ends at the centre (Hotulainen-Lappalainen
    # 2006 — required for the sarcomeric motor to rectify; "uniform" = all +z, e.g. filopodium/microvillus).
    graded = (fs.polarity == "mixed")
    fibers = []
    for f in range(n):
        zf = z[::-1] if (graded and f >= n // 2) else z
        fibers.append(np.column_stack([np.full(nb, xs[f]), np.full(nb, ys[f]), zf]))
    from ffn_sim.ff import units as U
    net = build_fiber_network(fibers, kappa=U.KAPPA_ACTIN)
    return net


def _build_lamellipodium_patch(spec: ArchitectureSpec, rng: np.random.Generator):
    """Arp2/3 DENDRITIC ±35° two-mode array on a flat patch (protrusion axis = +y, patch in x–y, z≈0).

    Mothers seeded at ±35° about +y; each spawns Arp2/3 daughters at the rest branch angle θ₀=70° toward
    the OTHER ∓35° mode (so every filament sits at ±35° two-mode AND every junction is at 70°). Returns
    (net, branch_triples (Nbr,3) [mother_after, branch_node, daughter_node1], anchor_pairs (Nbr,2)
    [daughter_node0, branch_node], fiber_angle (F,) signed orientation about +y)."""
    from ffn_sim.ff import units as U
    fs = spec.filament
    nb = max(3, int(round(fs.length_um / fs.seg_um)) + 1)
    seg, R = fs.seg_um, spec.R_um
    th0 = fs.branch_angle_rad
    sig = np.deg2rad(fs.branch_sigma_deg or 9.0)
    mode = np.deg2rad(fs.mode_axis_deg or 35.0)
    n = fs.n_filaments
    n_mothers = max(1, n // 2)

    def dir_of(phi):                                            # in-plane unit dir at angle φ from +y
        return np.array([np.sin(phi), np.cos(phi), 0.0])

    fibers, fiber_angle = [], []
    mothers = []                                               # (fiber_idx, phi)
    for m in range(n_mothers):
        s = 1.0 if rng.random() < 0.5 else -1.0
        phi = s * mode + sig * rng.standard_normal()
        base = np.array([rng.uniform(-R, R), rng.uniform(-R, R), 0.0])
        fibers.append(base[None, :] + np.arange(nb)[:, None] * seg * dir_of(phi)[None, :])
        fiber_angle.append(phi); mothers.append((m, s, phi))

    # daughters: branch off a random mother at a random interior bead, rotated by ∓θ₀ into the other mode
    dgt = []                                                   # (daughter_fiber_idx, mother_idx, branch_bead)
    for d in range(n - n_mothers):
        mi, s, mphi = mothers[rng.integers(n_mothers)]
        b = int(rng.integers(1, nb - 1))                       # interior branch bead on the mother
        dphi = mphi - s * th0 + sig * rng.standard_normal()    # → the other ∓35° mode; junction = θ₀
        branch_pos = fibers[mi][b]
        fibers.append(branch_pos[None, :] + np.arange(nb)[:, None] * seg * dir_of(dphi)[None, :])
        fiber_angle.append(dphi); dgt.append((n_mothers + d, mi, b))

    net = build_fiber_network(fibers, kappa=U.KAPPA_ACTIN)
    off = net.fiber_offsets
    triples, anchors = [], []
    for (di, mi, b) in dgt:
        ma = int(off[mi]) + min(b + 1, nb - 1)                 # mother node just past the branch
        bn = int(off[mi]) + b                                  # branch node (on mother)
        d0 = int(off[di])                                      # daughter base (co-located with bn)
        d1 = int(off[di]) + 1                                  # daughter first segment node
        triples.append([ma, bn, d1]); anchors.append([d0, bn])
    triples = np.array(triples, np.int64) if triples else np.zeros((0, 3), np.int64)
    anchors = np.array(anchors, np.int64) if anchors else np.zeros((0, 2), np.int64)
    return net, triples, anchors, np.array(fiber_angle)


def weave(spec: ArchitectureSpec, *, rng: np.random.Generator | None = None,
          alpha_fraction: float = 0.30) -> CrosslinkedCortex:
    """Build the woven network for ``spec`` → a CrosslinkedCortex (net + crosslinks + motors)."""
    if rng is None:
        rng = np.random.default_rng(0)
    fs = spec.filament
    n_xl = int(round(fs.n_filaments * spec.crosslinker.density_per_fil))
    n_myo = int(round(fs.n_filaments * spec.motor.density_per_fil)) if spec.motor.hand else 0

    branch_triples = np.zeros((0, 3), np.int64)
    anchors = np.zeros((0, 2), np.int64)
    xl_mask = myo_mask = None                                  # SF sarcomeric periodic node masks (set below)
    if spec.manifold == "sphere":
        # mirror gamma_floor.build_crosslinked_cortex (RNG order preserved → γ-floor parity)
        params = CortexParams(R_um=spec.R_um, n_filaments=fs.n_filaments,
                              beads_per_filament=int(round(fs.length_um / fs.seg_um)) + 1, seg_um=fs.seg_um)
        net, _ = build_cortex_network(params, rng=rng, n_filaments=fs.n_filaments)
        reach = mesoscale_reach(spec.R_um, fs.n_filaments)
    elif spec.manifold == "bundle":
        net = _build_bundle(spec, rng)
        spacing = fs.bundle_spacing_um if fs.bundle_spacing_um is not None else 0.008
        reach = 1.6 * spacing                                  # ~2× the bundle inter-filament spacing
        if fs.sarcomere_um:                                    # SF: periodic Z-body / anti-registered band planes
            z = net.pos[:, 2]; p = fs.sarcomere_um; w = 0.5 * fs.seg_um
            xl_mask = np.abs((z + p / 2) % p - p / 2) < w      # α-actinin Z-bodies @ k·period
            myo_mask = np.abs(z % p - p / 2) < w               # NMIIA bands @ (k+½)·period (anti-registered)
    elif spec.manifold == "patch":                            # lamellipodium dendritic array (Arp2/3)
        net, branch_triples, anchors, _ = _build_lamellipodium_patch(spec, rng)
        reach = 0.6                                            # lamellipodial mesh ξ ~0.6 µm (Sakamoto 2024)
    else:
        raise ValueError(f"manifold {spec.manifold!r} not supported (sphere|bundle|patch)")

    # crosslinkers: KDTree near cross-fiber pairs, filtered by bind mode (+ SF Z-body periodic mask)
    xl_pairs = _cross_fiber_pairs(net, reach, n_xl, rng, node_mask=xl_mask)
    xl_pairs = _filter_bind_mode(net, xl_pairs, spec.crosslinker.bind_mode)
    nxl = xl_pairs.shape[0]
    if spec.manifold == "sphere":                              # cortex α-actinin/filamin split (parity)
        is_alpha = rng.random(nxl) < alpha_fraction
        xl_k = np.where(is_alpha, ALPHA_ACTININ.link_k, FILAMIN.link_k)
    else:                                                      # bundle/patch: single crosslinker (spec.hand)
        xl_k = np.full(nxl, spec.crosslinker.hand.link_k)
    xl_i, xl_j = (xl_pairs[:, 0], xl_pairs[:, 1]) if nxl else (np.zeros(0, np.int64), np.zeros(0, np.int64))
    xl_rest = np.linalg.norm(net.pos[xl_j] - net.pos[xl_i], axis=1) if nxl else np.zeros(0)

    if anchors.shape[0]:                                       # Arp2/3 branch ANCHORS (daughter base ↔ branch
        a_i, a_j = anchors[:, 0], anchors[:, 1]               # node): stiff, force-free at the branch geometry
        a_rest = np.linalg.norm(net.pos[a_j] - net.pos[a_i], axis=1)
        xl_i = np.concatenate([xl_i, a_i]); xl_j = np.concatenate([xl_j, a_j])
        xl_k = np.concatenate([xl_k, np.full(anchors.shape[0], ALPHA_ACTININ.link_k)])
        xl_rest = np.concatenate([xl_rest, a_rest])

    # motors (SF: anti-registered band-centre periodic mask → sarcomeric placement)
    if n_myo:
        used = {(int(a), int(b)) for a, b in xl_pairs}
        myo_pairs = _cross_fiber_pairs(net, reach, n_myo, rng, exclude=used, node_mask=myo_mask)
    else:
        myo_pairs = np.zeros((0, 2), np.int64)

    r0_mean = float(np.linalg.norm(net.pos - net.pos.mean(axis=0), axis=1).mean())
    return CrosslinkedCortex(
        net=net, xl_i=xl_i, xl_j=xl_j, xl_k=xl_k, xl_rest=xl_rest,
        myo_i=myo_pairs[:, 0], myo_j=myo_pairs[:, 1], R_um=spec.R_um, R0_mean=r0_mean,
        branch_triples=branch_triples)
