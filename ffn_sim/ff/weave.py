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
    spacing_um = 0.008                                          # ~8 nm inter-filament (fascin geometry; PI-gated)
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
    fibers = [np.column_stack([np.full(nb, xs[f]), np.full(nb, ys[f]), z]) for f in range(n)]
    from ffn_sim.ff import units as U
    net = build_fiber_network(fibers, kappa=U.KAPPA_ACTIN)
    return net


def weave(spec: ArchitectureSpec, *, rng: np.random.Generator | None = None,
          alpha_fraction: float = 0.30) -> CrosslinkedCortex:
    """Build the woven network for ``spec`` → a CrosslinkedCortex (net + crosslinks + motors)."""
    if rng is None:
        rng = np.random.default_rng(0)
    fs = spec.filament
    n_xl = int(round(fs.n_filaments * spec.crosslinker.density_per_fil))
    n_myo = int(round(fs.n_filaments * spec.motor.density_per_fil)) if spec.motor.hand else 0

    if spec.manifold == "sphere":
        # mirror gamma_floor.build_crosslinked_cortex (RNG order preserved → γ-floor parity)
        params = CortexParams(R_um=spec.R_um, n_filaments=fs.n_filaments,
                              beads_per_filament=int(round(fs.length_um / fs.seg_um)) + 1, seg_um=fs.seg_um)
        net, _ = build_cortex_network(params, rng=rng, n_filaments=fs.n_filaments)
        reach = mesoscale_reach(spec.R_um, fs.n_filaments)
    elif spec.manifold == "bundle":
        net = _build_bundle(spec, rng)
        reach = 1.6 * 0.008                                    # ~2× the 8 nm bundle spacing
    else:
        raise ValueError(f"manifold {spec.manifold!r} not in increment-1 scope (sphere|bundle)")

    # crosslinkers: KDTree near cross-fiber pairs, filtered by bind mode
    xl_pairs = _cross_fiber_pairs(net, reach, n_xl, rng)
    xl_pairs = _filter_bind_mode(net, xl_pairs, spec.crosslinker.bind_mode)
    nxl = xl_pairs.shape[0]
    if spec.manifold == "sphere":                              # cortex α-actinin/filamin split (parity)
        is_alpha = rng.random(nxl) < alpha_fraction
        xl_k = np.where(is_alpha, ALPHA_ACTININ.link_k, FILAMIN.link_k)
    else:                                                      # bundle: single crosslinker (spec.hand)
        xl_k = np.full(nxl, spec.crosslinker.hand.link_k)
    xl_rest = (np.linalg.norm(net.pos[xl_pairs[:, 1]] - net.pos[xl_pairs[:, 0]], axis=1)
               if nxl else np.zeros(0))

    # motors
    if n_myo:
        used = {(int(a), int(b)) for a, b in xl_pairs}
        myo_pairs = _cross_fiber_pairs(net, reach, n_myo, rng, exclude=used)
    else:
        myo_pairs = np.zeros((0, 2), np.int64)

    r0_mean = float(np.linalg.norm(net.pos - net.pos.mean(axis=0), axis=1).mean())
    return CrosslinkedCortex(
        net=net, xl_i=xl_pairs[:, 0], xl_j=xl_pairs[:, 1], xl_k=xl_k, xl_rest=xl_rest,
        myo_i=myo_pairs[:, 0], myo_j=myo_pairs[:, 1], R_um=spec.R_um, R0_mean=r0_mean)
