"""Connected contractile fiber-network — buckling-enabled active-stress measurement (Stage 6h).

The decisive test from FF_STAGE6D §3.0c: does a buckling-capable, percolating, contractile cortex
network develop band-level active stress, or does the floor survive? This module builds a connected
inextensible contractile network (the robust solver the prior attempt lacked) and measures the
macroscopic contractile stress σ as a function of the per-link contractile force and the network
connectivity z.

Construction (all geometry-fixed or lit-anchored — connectivity is a CONTROLLED VARIABLE, swept, not
tuned):
  * a 2D triangular lattice (interior z=6) of vertices, each lattice EDGE a 3-node fiber (2 segments
    → can buckle); diluting edges (keep-fraction) sweeps the mean coordination z down through the
    sub-isostatic regime (z<4 in 2D central-force) where Ronceray/Broedersz/Lenz 2016 amplification,
    if it operates, would appear;
  * fibers welded at shared vertices by stiff crosslink springs (percolating network, giant
    component only);
  * a contractile motor on each fiber (constant force f_act pulling its endpoints together — the
    active prestress; swept as the controlled variable);
  * boundary vertices pinned; relaxed by explicit projected-free descent + periodic reshape
    (robust — no singular projector, no M-SHAKE; the FF engine buckles, `test_buckling.py`).

σ is read as the contractile line tension = boundary reaction / perimeter [pN/µm] (1 pN/µm = 1e-3
mN/m; band = 350–650 pN/µm). FINDING (sweeps, 2026-06-30): σ ≈ 10–90 pN/µm at f_act 1–100 pN, ~30–100×
UNDER band, and **flat in z across 2.98–5.29** — connectivity/buckling does NOT lift the floor in the
INEXTENSIBLE network. (Ronceray's ~7× amplification is an EXTENSIBLE-network effect; FF/Cytosim
inextensibility may structurally preclude it — a finding surfaced to PI, not a tuned knob.)
"""

from __future__ import annotations

import collections
from dataclasses import dataclass

import numpy as np

from aleph.laws import units as U
from aleph.laws.constraints import reshape
from aleph.laws.fiber_network import build_fiber_network
from aleph.laws.forces_warp import make_bending_force_fn


@dataclass(slots=True)
class ContractilePatch:
    """A built connected contractile network patch + its solver handles."""

    net: object
    xl: np.ndarray            # (Nxl, 2) crosslink weld-spring node pairs (rest 0)
    ep: list                  # [(i0, i1)] contractile-motor endpoints per fiber
    pin: np.ndarray           # (N,) bool pinned (boundary) nodes
    coord: np.ndarray         # lattice vertex coords
    z_mean: float             # mean coordination of the giant component
    bbox: tuple               # (xmin, xmax, ymin, ymax)


def build_patch(*, nx: int = 11, ny: int = 11, a: float = 0.3, kappa: float = U.KAPPA_ACTIN,
                keep_frac: float = 1.0, seed_z: float = 0.08,
                rng: np.random.Generator | None = None) -> ContractilePatch | None:
    """Build a (diluted) triangular-lattice contractile network; None if no percolating component."""
    if rng is None:
        rng = np.random.default_rng(0)
    vert = {}
    coord = []
    idx = 0
    for j in range(ny):
        for i in range(nx):
            vert[(i, j)] = idx
            coord.append(((i + 0.5 * (j % 2)) * a, j * a * np.sqrt(3) / 2, 0.0))
            idx += 1
    coord = np.array(coord)
    edges = set()
    for j in range(ny):
        for i in range(nx):
            for di, dj in [(1, 0), (0, 1), (-1, 1) if j % 2 == 0 else (1, 1)]:
                if (i + di, j + dj) in vert:
                    A, B = vert[(i, j)], vert[(i + di, j + dj)]
                    edges.add((min(A, B), max(A, B)))
    keep = [e for e in sorted(edges) if rng.random() < keep_frac]
    adj = collections.defaultdict(list)
    for A, B in keep:
        adj[A].append(B); adj[B].append(A)
    seen = set(); comps = []
    for s in list(adj):
        if s in seen:
            continue
        stack = [s]; comp = []
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            seen.add(u); comp.append(u); stack += [w for w in adj[u] if w not in seen]
        comps.append(comp)
    if not comps:
        return None
    giant = set(max(comps, key=len))
    keep = [(A, B) for A, B in keep if A in giant and B in giant]
    if len(keep) < 5:
        return None
    deg = collections.Counter()
    for A, B in keep:
        deg[A] += 1; deg[B] += 1
    fibers = []
    fev = []
    for (A, B) in keep:
        p0, p1 = coord[A], coord[B]
        mid = 0.5 * (p0 + p1) + np.array([0, 0, seed_z * a * rng.standard_normal()])
        fibers.append(np.array([p0, mid, p1])); fev.append((A, B))
    net = build_fiber_network(fibers, kappa=kappa)
    off = net.fiber_offsets
    vn = {}
    for f, (A, B) in enumerate(fev):
        vn.setdefault(A, []).append(int(off[f])); vn.setdefault(B, []).append(int(off[f + 1]) - 1)
    xl = [(nodes[0], nodes[k]) for nodes in vn.values() for k in range(1, len(nodes))]
    xl = np.array(xl) if xl else np.zeros((0, 2), int)
    ep = [(int(off[f]), int(off[f + 1]) - 1) for f in range(net.n_fibers)]
    xy = coord[:, :2]
    xmin, xmax, ymin, ymax = xy[:, 0].min(), xy[:, 0].max(), xy[:, 1].min(), xy[:, 1].max()
    pin = np.zeros(net.n_nodes, bool)
    for v, nodes in vn.items():
        x, y = coord[v, 0], coord[v, 1]
        if x < xmin + 1e-6 or x > xmax - 1e-6 or y < ymin + 1e-6 or y > ymax - 1e-6:
            for nd in nodes:
                pin[nd] = True
    if pin.sum() < 3:
        return None
    z_mean = float(np.mean([deg[v] for v in giant if deg[v] > 0]))
    return ContractilePatch(net=net, xl=xl, ep=ep, pin=pin, coord=coord, z_mean=z_mean,
                            bbox=(xmin, xmax, ymin, ymax))


def _forces(patch: ContractilePatch, x: np.ndarray, bfn, f_act: float, k_xl: float) -> np.ndarray:
    N = patch.net.n_nodes
    F = bfn(x.reshape(-1)).reshape(N, 3).copy()
    if patch.xl.shape[0]:
        d = x[patch.xl[:, 1]] - x[patch.xl[:, 0]]
        f0 = k_xl * d
        np.add.at(F, patch.xl[:, 0], f0); np.add.at(F, patch.xl[:, 1], -f0)
    for (i0, i1) in patch.ep:
        d = x[i1] - x[i0]
        L = float(np.linalg.norm(d))
        if L > 1e-9:
            u = d / L
            F[i0] += f_act * u; F[i1] -= f_act * u
    return F


def measure_sigma(patch: ContractilePatch, f_act: float, *, a: float = 0.3,
                  kappa: float = U.KAPPA_ACTIN, k_xl: float = 20.0, n_steps: int = 9000) -> dict | None:
    """Relax the contractile patch and return the macroscopic contractile stress σ [pN/µm] + diagnostics.

    σ = boundary reaction · r̂ / perimeter (the contractile line tension the patch develops). Also
    returns the raw dipole-sum baseline σ_dipole and the out-of-plane buckling amplitude.
    """
    net = patch.net
    N = net.n_nodes
    pos = net.pos.copy()
    pinpos = pos[patch.pin].copy()
    dt_mu = 0.2 / max(16 * kappa / a**3, k_xl)
    bfn = make_bending_force_fn(net)
    x = pos.copy()
    for it in range(n_steps):
        net.pos = x
        F = _forces(patch, x, bfn, f_act, k_xl)
        F[patch.pin] = 0.0
        if not np.isfinite(F).all():
            return None
        x = x + dt_mu * F
        if (it + 1) % 100 == 0:
            net.pos = x
            xr = reshape(net, n_iter=1)
            xr[patch.pin] = pinpos
            x = xr
    net.pos = x
    Ff = _forces(patch, x, bfn, f_act, k_xl)
    react = -Ff[patch.pin]
    cen = x[~patch.pin].mean(axis=0)
    rxy = (x[patch.pin] - cen)[:, :2]
    rn = np.linalg.norm(rxy, axis=1) + 1e-12
    xmin, xmax, ymin, ymax = patch.bbox
    perim = 2 * ((xmax - xmin) + (ymax - ymin))
    sigma = float((react[:, :2] * rxy / rn[:, None]).sum(1).sum() / perim)   # +contractile
    Ldip = np.array([np.linalg.norm(x[i1] - x[i0]) for i0, i1 in patch.ep])
    sigma_dipole = float((f_act * Ldip).sum() / ((xmax - xmin) * (ymax - ymin)) / 2)
    return {"sigma_pN_um": sigma, "sigma_dipole_pN_um": sigma_dipole, "z_mean": patch.z_mean,
            "zbow_um": float(np.max(np.abs(x[:, 2]))), "f_act": f_act}
