"""ECM as a 3D Mikado fiber network (collagen-I) — the extracellular matrix a protrusion pushes into.

PI 2026-07-02 ("ecm도 mikado 네트워크로 구성해서"): the ECM is not a smooth wall or fixed anchor points
(that is the FA-clutch substrate, ``fa_anchor``) — it is a disordered network of stiff collagen fibers. The
**Mikado model** (Wilhelm & Frey 2003, PRL 91:108103; Head-Levine-MacKintosh 2003, PRE 68:061907) is the
canonical random-fiber network: straight rods placed at random positions + orientations, cross-linked where
they (nearly) intersect. This is the 3D collagen-gel generalisation — random rods in a slab, cross-linked at
near-contacts, with the far boundary pinned so the network is embedded in bulk ECM and RESISTS a protrusion.

Mechanistic, reuses the FF force kernels unchanged: each collagen fiber is a Cytosim discrete fiber
(``fiber_network`` — segments + bending triples), cross-links are Hookean ``link_spring`` bonds at
near-contacts, boundary anchoring is a pinned-node BC. The protrusion↔ECM interaction (soft excluded-volume
contact) is applied by the driver via ``network_warp.soft_contact_kernel`` — the ECM then supplies the
protrusion's load EMERGENTLY (it pushes fibers aside / stalls), not as an imposed number.

Constants (units: µm, pN, s — FF convention):
- collagen fiber persistence length ``LP_COLLAGEN_UM`` → κ = kBT·Lp. Reconstituted collagen-I fibrils span a
  wide Lp (thin fibril ~tens of µm; thick bundles → mm); we use a mid thin-fibril value, **PI-gated** (no single
  MCF7-matrix datum). κ_collagen ≫ κ_actin, so the ECM is much stiffer than the cortex — the physical point.
- mesh/pore size ξ: an EMERGENT geometric output of (fiber length density) here, reported by the builder;
  reconstituted collagen-I gels are ξ ~ 1–5 µm (e.g. Yang & Kaufman 2009, Biophys J 96:1566; 2.5 mg/mL ≈ 2–3 µm)
  — a target the builder is set to hit, not a tuned mechanism.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from ffn_sim.ff import units as U
from ffn_sim.ff.fiber_network import FiberNetwork, build_fiber_network

# Collagen-I fiber bending — PI-gated (documented range, not a fitted value). Thin reconstituted collagen
# fibril persistence length; κ = kBT·Lp (same law as KAPPA_ACTIN). ~20 µm ⇒ κ ≈ 0.086 pN·µm² … but collagen
# BUNDLES are far stiffer; this is the compliant end. Flagged for PI KB-registration (collagen Lp source).
LP_COLLAGEN_UM = 20.0
KAPPA_COLLAGEN = U.KBT * LP_COLLAGEN_UM            # pN·µm²

# Collagen cross-link stiffness — entanglement / covalent cross-link; PI-gated. Uses a stiff Hookean bond in
# the range of the actin crosslinkers (network response is dominated by geometry + pinning, not this value).
K_XL_COLLAGEN = 100.0                              # pN/µm  (PI-gated; documented, not tuned to an outcome)


@dataclass(slots=True)
class MikadoNetwork:
    """A 3D Mikado ECM network: discrete collagen fibers + near-contact cross-links + a pinned boundary."""

    net: FiberNetwork
    xl_i: np.ndarray            # (X,) crosslink node i
    xl_j: np.ndarray            # (X,) crosslink node j (different fiber, within contact radius)
    xl_k: np.ndarray            # (X,) crosslink stiffness [pN/µm]
    xl_rest: np.ndarray         # (X,) crosslink rest length [µm]
    pinned: np.ndarray          # (N,) bool — boundary-anchored nodes (Dirichlet BC; embedded in bulk ECM)
    box_lo: np.ndarray          # (3,) µm
    box_hi: np.ndarray          # (3,) µm
    mesh_size_um: float         # emergent pore size ξ (median inter-crosslink spacing along fibers)


def _random_rod(box_lo, box_hi, length, rng):
    """A random straight rod: uniform centre in box, isotropic orientation; endpoints clipped to the box."""
    c = box_lo + rng.random(3) * (box_hi - box_lo)
    v = rng.standard_normal(3)
    v /= np.linalg.norm(v) + 1e-12
    a, b = c - 0.5 * length * v, c + 0.5 * length * v
    # clip both ends into the box (keep the rod inside; may shorten it near faces)
    a = np.clip(a, box_lo, box_hi)
    b = np.clip(b, box_lo, box_hi)
    return a, b


def build_mikado_network(box_lo, box_hi, *, n_fibers: int, fiber_len_um: float = 5.0, seg_um: float = 0.5,
                         kappa: float = KAPPA_COLLAGEN, xl_contact_um: float = 0.25,
                         xl_k: float = K_XL_COLLAGEN, pin_face: str = "x_lo", pin_margin_um: float = 0.5,
                         rng: np.random.Generator | None = None) -> MikadoNetwork:
    """Build a 3D Mikado collagen ECM network in ``[box_lo, box_hi]`` (µm).

    Args:
        box_lo, box_hi: (3,) box corners [µm].
        n_fibers: number of collagen rods (sets the length density → mesh size ξ).
        fiber_len_um: rod contour length L [µm].
        seg_um: discretisation segment length ℓ₀ [µm].
        kappa: collagen fiber bending rigidity [pN·µm²] (κ = kBT·Lp; PI-gated).
        xl_contact_um: two nodes on DIFFERENT fibers closer than this are cross-linked (the Mikado
            near-intersection rule; ~ collagen fiber diameter scale).
        xl_k: cross-link Hookean stiffness [pN/µm].
        pin_face: which box face is Dirichlet-pinned (embedded in bulk ECM), e.g. ``"x_lo"``; nodes within
            ``pin_margin_um`` of that face are fixed. Use ``"none"`` for a free network.
        pin_margin_um: pin-band thickness [µm].
        rng: numpy Generator.

    Returns:
        A :class:`MikadoNetwork` (fiber network + near-contact cross-links + pinned BC mask + mesh size).
    """
    if rng is None:
        rng = np.random.default_rng(0)
    box_lo = np.asarray(box_lo, np.float64)
    box_hi = np.asarray(box_hi, np.float64)
    nb = max(2, int(round(fiber_len_um / seg_um)) + 1)

    fibers = []
    for _ in range(n_fibers):
        a, b = _random_rod(box_lo, box_hi, fiber_len_um, rng)
        if np.linalg.norm(b - a) < seg_um:                 # degenerate (clipped to ~point) → resample once
            a, b = _random_rod(box_lo, box_hi, fiber_len_um, rng)
        t = np.linspace(0.0, 1.0, nb)[:, None]
        fibers.append(a[None, :] * (1 - t) + b[None, :] * t)
    net = build_fiber_network(fibers, kappa=kappa)

    # cross-links at near-contacts between DIFFERENT fibers (the Mikado intersection rule, 3D near-contact)
    node_fiber = np.concatenate([np.full(net.fiber_offsets[f + 1] - net.fiber_offsets[f], f)
                                 for f in range(net.n_fibers)])
    tree = cKDTree(net.pos)
    pairs = tree.query_pairs(r=xl_contact_um, output_type="ndarray")
    if pairs.shape[0]:
        diff_fiber = node_fiber[pairs[:, 0]] != node_fiber[pairs[:, 1]]
        pairs = pairs[diff_fiber]
    # keep at most one crosslink per fiber-pair (nearest contact) so we don't over-stitch
    if pairs.shape[0]:
        keyfp = node_fiber[pairs[:, 0]].astype(np.int64) * net.n_fibers + node_fiber[pairs[:, 1]]
        d = np.linalg.norm(net.pos[pairs[:, 1]] - net.pos[pairs[:, 0]], axis=1)
        order = np.argsort(d)
        seen, keep = set(), []
        for idx in order:
            k = int(keyfp[idx])
            if k in seen:
                continue
            seen.add(k); keep.append(idx)
        pairs = pairs[np.array(keep)]
    xl_i = pairs[:, 0] if pairs.shape[0] else np.zeros(0, np.int64)
    xl_j = pairs[:, 1] if pairs.shape[0] else np.zeros(0, np.int64)
    xl_rest = np.linalg.norm(net.pos[xl_j] - net.pos[xl_i], axis=1) if pairs.shape[0] else np.zeros(0)
    xl_k_arr = np.full(xl_i.shape[0], float(xl_k))

    # boundary pin (Dirichlet BC): nodes within pin_margin of the named face → embedded in bulk ECM
    pinned = np.zeros(net.n_nodes, bool)
    if pin_face != "none":
        ax = {"x": 0, "y": 1, "z": 2}[pin_face[0]]
        if pin_face.endswith("lo"):
            pinned = net.pos[:, ax] < box_lo[ax] + pin_margin_um
        else:
            pinned = net.pos[:, ax] > box_hi[ax] - pin_margin_um

    # mesh size ξ = median inter-crosslink spacing along fibers (practical pore scale)
    if xl_i.shape[0]:
        xl_nodes = np.concatenate([xl_i, xl_j])
        # spacing between crosslink nodes that lie on the same fiber, sorted by arclength (node index proxy)
        spac = []
        for f in range(net.n_fibers):
            a0, b0 = net.fiber_offsets[f], net.fiber_offsets[f + 1]
            on = np.sort(xl_nodes[(xl_nodes >= a0) & (xl_nodes < b0)])
            if on.size >= 2:
                spac.extend(np.diff(on) * seg_um)
        mesh = float(np.median(spac)) if spac else fiber_len_um
    else:
        mesh = fiber_len_um

    return MikadoNetwork(net=net, xl_i=xl_i, xl_j=xl_j, xl_k=xl_k_arr, xl_rest=xl_rest,
                         pinned=pinned, box_lo=box_lo, box_hi=box_hi, mesh_size_um=mesh)
