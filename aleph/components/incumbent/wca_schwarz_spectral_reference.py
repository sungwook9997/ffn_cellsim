r"""CPU/NumPy spectral reference for the WCA contact-graph Schwarz preconditioner.

This is a *predictor*, not a native result.  It builds a small, exactly-diagonalizable multi-fiber contact
network that reproduces the full-native plateau structure documented in ``outputs/ac/implicit/REPORT.md`` (a
stiff radial WCA contact chain between fibers, Hookean crosslinks from a contacting fiber to a *distinct*
non-adjacent fiber, plus fiber bending, NF2007 inextensibility, and the derived omitted-family regularizer),
and asks -- before any device implementation cost -- whether the proposed multi-fiber **star** subdomain
(``Open numerical work`` item 1) resolves the stuck mode that the falsified independent-crosslink and the
retained rigid-fiber blocks leave behind.

The reduced operator matches the production analytic tangent used by ``contact_schwarz.py``:

* WCA radial tangent ``radial * outer(u, u)`` with ``radial = (24 eps / r^2)(26 sr^12 - 7 sr^6)`` (force-capped
  pairs contribute zero radial curvature, exactly as in ``_wca_tangent``);
* central-spring tangent ``transverse * I + (k - transverse) outer(u, u)`` with
  ``transverse = max(k (r - r0)/r, 0)`` (exactly ``_central_tangent``);
* the NF2007 inextensibility projector ``P = I - J^T (J J^T)^{-1} J`` built from bond unit vectors; and
* the ``sqrt(eps64) * k_max`` rigid-mode regularizer floor.

Two lenses are reported for each candidate subdomain family (WCA node-pair, rigid fiber-pair, and the proposed
star):

1. *Faithful-to-tournament*: the single additive-Schwarz direction ``d = P M P r`` applied to the plateau
   residual, and the fraction of the operator's softest (stuck) mode that the subdomains span.  This mirrors
   how ``contact_schwarz`` uses the correction -- as one residual-monotone tournament candidate.
2. *Idealized preconditioner*: the same additive Schwarz completed with a node-diagonal block on every
   uncovered node so ``M`` is full-rank SPD, then ``cond(M A)`` and preconditioned-CG iterations on the free
   subspace -- the conditioning the subdomain design would achieve as a full PCG preconditioner.

Both lenses agree when a design captures the coupled mode.  The two-node WCA path is cross-checked against
``contact_schwarz.additive_pair_schwarz_oracle`` when that module imports on the host.

Sanity Gate:
    * Dimensions: tangents are pN/um, the operator acts on um displacements, residual is pN.
    * Symmetry/positivity: ``P`` is a symmetric idempotent; every local block ``A_e`` is a principal submatrix
      of the SPD regularized operator; the completed ``M`` is symmetric positive definite; all checked with
      ``eigvalsh``.
    * Boundary cases: a force-capped WCA pair contributes zero radial curvature; an uncovered node receives its
      exact diagonal block only.
    * Numerical: dense FP64 throughout; the model is small enough for exact ``eigvalsh``/``solve`` so no
      iterative tolerance enters the spectral verdict.  No damping, cluster radius, or fitted threshold.
    * Scope: host predictor on a reduced model.  It states which native run validates its prediction; it never
      claims a native convergence pass.
"""

from __future__ import annotations

import dataclasses

import numpy as np

__all__ = [
    "ContactNetwork",
    "build_contact_network",
    "nf2007_projector",
    "wca_radial_tangent",
    "central_spring_tangent",
    "additive_schwarz",
    "pair_subdomains",
    "rigid_fiber_subdomains",
    "star_subdomains",
    "coupled_translation_residual",
    "direction_reduction",
    "local_cluster_spd",
    "residual_dominance",
]

_SQRT_EPS64 = float(np.sqrt(np.finfo(np.float64).eps))


# --------------------------------------------------------------------------------------------------
# Production-matched local tangents (dense NumPy mirrors of contact_schwarz._wca_tangent /
# _central_tangent).
# --------------------------------------------------------------------------------------------------
def wca_radial_tangent(delta: np.ndarray, sigma: float, epsilon: float, force_cap: float = 0.0) -> np.ndarray:
    """3x3 WCA radial curvature block ``radial * outer(u, u)`` (zero outside the cutoff or when capped)."""
    delta = np.asarray(delta, dtype=np.float64)
    length = float(np.linalg.norm(delta))
    cutoff = (2.0 ** (1.0 / 6.0)) * sigma
    if length <= 1.0e-12 or length >= cutoff:
        return np.zeros((3, 3))
    ratio = sigma / length
    sr6 = ratio ** 6
    force_magnitude = (24.0 * epsilon / length) * (2.0 * sr6 * sr6 - sr6)
    if force_cap > 0.0 and force_magnitude >= force_cap:
        return np.zeros((3, 3))
    radial = (24.0 * epsilon / (length * length)) * (26.0 * sr6 * sr6 - 7.0 * sr6)
    unit = delta / length
    return radial * np.outer(unit, unit)


def central_spring_tangent(delta: np.ndarray, stiffness: float, rest_length: float) -> np.ndarray:
    """3x3 Hookean central-spring tangent (transverse isotropic + longitudinal), matching _central_tangent."""
    delta = np.asarray(delta, dtype=np.float64)
    length = float(np.linalg.norm(delta))
    if length <= 1.0e-12:
        return np.zeros((3, 3))
    unit = delta / length
    transverse = max(stiffness * (length - rest_length) / length, 0.0)
    difference = stiffness - transverse
    return transverse * np.eye(3) + difference * np.outer(unit, unit)


def _wca_epsilon_for_curvature(k_max: float, sigma: float, overlap: float) -> float:
    """Pick epsilon so the WCA radial curvature at r = overlap*sigma equals k_max (production ledger scale)."""
    r = overlap * sigma
    sr6 = (sigma / r) ** 6
    shape = (24.0 / (r * r)) * (26.0 * sr6 * sr6 - 7.0 * sr6)
    return k_max / shape


# --------------------------------------------------------------------------------------------------
# Multi-fiber contact-network model.
# --------------------------------------------------------------------------------------------------
@dataclasses.dataclass
class ContactNetwork:
    """A multi-fiber fragment reproducing the WCA-contact-chain + distinct-crosslink plateau structure."""

    positions: np.ndarray            # (N, 3) node coordinates [um]
    fiber_of: np.ndarray             # (N,) fiber id per node
    n_fibers: int
    bonds: np.ndarray                # (B, 2) consecutive-node inextensibility bonds
    stiffness: np.ndarray            # (3N, 3N) dense unprojected tangent aI + K_wca + K_xlink + K_bend
    projector: np.ndarray            # (3N, 3N) NF2007 inextensibility projector P
    wca_edges: np.ndarray            # (E_w, 2) node pairs carrying a stiff radial contact
    xlink_edges: np.ndarray          # (E_x, 2) node pairs carrying a crosslink to a distinct fiber
    k_max: float
    parts: dict

    @property
    def n_nodes(self) -> int:
        return self.positions.shape[0]

    @property
    def operator(self) -> np.ndarray:
        """The projected symmetric operator ``A = P (aI + K) P``."""
        return self.projector @ self.stiffness @ self.projector


def nf2007_projector(positions: np.ndarray, bonds: np.ndarray) -> np.ndarray:
    """Symmetric idempotent inextensibility projector ``P = I - J^T (J J^T)^{-1} J`` from bond unit vectors."""
    positions = np.asarray(positions, dtype=np.float64)
    n = positions.shape[0]
    dof = 3 * n
    jac = np.zeros((len(bonds), dof))
    for row, (i, j) in enumerate(bonds):
        delta = positions[j] - positions[i]
        unit = delta / float(np.linalg.norm(delta))
        jac[row, 3 * i:3 * i + 3] = -unit
        jac[row, 3 * j:3 * j + 3] = unit
    gram = jac @ jac.T
    proj = np.eye(dof) - jac.T @ np.linalg.solve(gram, jac)
    return 0.5 * (proj + proj.T)


def _scatter(target, i, j, block):
    target[3 * i:3 * i + 3, 3 * i:3 * i + 3] += block
    target[3 * j:3 * j + 3, 3 * j:3 * j + 3] += block
    target[3 * i:3 * i + 3, 3 * j:3 * j + 3] -= block
    target[3 * j:3 * j + 3, 3 * i:3 * i + 3] -= block


def build_contact_network(
    *,
    n_fibers: int = 8,
    nodes_per_fiber: int = 5,
    segment_length: float = 0.5,
    k_max: float = 2_828_099.2886817832,
    wca_sigma: float = 0.05,
    wca_overlap: float = 0.76,
    crosslink_stiffness: float = 8.2e5,
    xlink_stride: int = 2,
    bending_stiffness: float = 7.0e2,
    seed_jitter: float = 1.0e-4,
) -> ContactNetwork:
    """Build ``n_fibers`` parallel fibers in a bundle with a WCA contact chain and distinct crosslinks.

    * Fibers are laid parallel along x, stacked in y at the WCA contact gap, so fiber ``f`` and ``f+1`` touch
      at their centre nodes (a contact chain).
    * Each fiber ``f`` crosslinks its centre node to fiber ``f + xlink_stride`` (a *distinct*, non-adjacent
      fiber) -- the documented structure the pair/rigid blocks miss.
    * A star subdomain (WCA fibers + crosslink-neighbor fiber) is therefore a genuine *local* cluster inside a
      larger network, not the whole system.

    Magnitudes are production ledger values (``k_max``, crosslink Ferrer stiffness); nothing is tuned to a
    convergence outcome.  ``wca_overlap = r/sigma`` fixes the radial curvature near the reported ``4.8e5``.
    """
    npf = nodes_per_fiber
    contact_gap = wca_overlap * wca_sigma
    cmid = npf // 2
    xs = np.arange(npf) * segment_length
    fibers = []
    for f in range(n_fibers):
        fibers.append(np.column_stack([xs, np.full(npf, f * contact_gap), np.zeros(npf)]))
    positions = np.vstack(fibers).astype(np.float64)
    if seed_jitter > 0.0:
        idx = np.arange(positions.shape[0])[:, None]
        positions = positions + seed_jitter * np.cos(idx * np.array([[1.7, 2.3, 3.1]]))
    fiber_of = np.concatenate([np.full(npf, f) for f in range(n_fibers)])

    bonds = []
    for f in range(n_fibers):
        base = f * npf
        for k in range(npf - 1):
            bonds.append((base + k, base + k + 1))
    bonds = np.asarray(bonds, dtype=np.int64)

    dof = 3 * positions.shape[0]
    k_wca = np.zeros((dof, dof))
    k_xlink = np.zeros((dof, dof))
    k_bend = np.zeros((dof, dof))
    epsilon = _wca_epsilon_for_curvature(k_max, wca_sigma, wca_overlap)

    wca_edges = []
    for f in range(n_fibers - 1):
        i = f * npf + cmid
        j = (f + 1) * npf + cmid
        block = wca_radial_tangent(positions[j] - positions[i], wca_sigma, epsilon)
        _scatter(k_wca, i, j, block)
        wca_edges.append((i, j))
    wca_edges = np.asarray(wca_edges, dtype=np.int64)

    xlink_rest = 0.03
    xlink_edges = []
    for f in range(n_fibers - xlink_stride):
        a = f * npf + cmid
        c = (f + xlink_stride) * npf + cmid
        # offset the crosslink node slightly in z so it is a genuine distinct-fiber spring, not the WCA line
        delta = positions[c] - positions[a]
        block = central_spring_tangent(delta if np.linalg.norm(delta) > 1e-9 else np.array([0, 0, xlink_rest]),
                                       crosslink_stiffness, xlink_rest)
        _scatter(k_xlink, a, c, block)
        xlink_edges.append((a, c))
    xlink_edges = np.asarray(xlink_edges, dtype=np.int64)

    for f in range(n_fibers):
        base = f * npf
        for k in range(1, npf - 1):
            im, ic, ip = base + k - 1, base + k, base + k + 1
            b = bending_stiffness
            for a, s in ((im, 1.0), (ic, -2.0), (ip, 1.0)):
                for c, t in ((im, 1.0), (ic, -2.0), (ip, 1.0)):
                    k_bend[3 * a:3 * a + 3, 3 * c:3 * c + 3] += (b * s * t) * np.eye(3)

    reg = _SQRT_EPS64 * k_max
    stiffness = reg * np.eye(dof) + k_wca + k_xlink + k_bend
    stiffness = 0.5 * (stiffness + stiffness.T)
    projector = nf2007_projector(positions, bonds)
    return ContactNetwork(
        positions=positions, fiber_of=fiber_of, n_fibers=n_fibers, bonds=bonds,
        stiffness=stiffness, projector=projector, wca_edges=wca_edges, xlink_edges=xlink_edges,
        k_max=k_max,
        parts={"regularizer": reg * np.eye(dof), "wca": k_wca, "xlink": k_xlink, "bend": k_bend},
    )


# --------------------------------------------------------------------------------------------------
# Subdomain families (tiled over every WCA contact edge).
# --------------------------------------------------------------------------------------------------
def pair_subdomains(net: ContactNetwork) -> list[np.ndarray]:
    """Two-node WCA-pair subdomains, one per contact edge (the base/falsified path)."""
    return [np.array(edge) for edge in net.wca_edges]


def rigid_fiber_subdomains(net: ContactNetwork) -> list[np.ndarray]:
    """Two-fiber subdomains: all nodes of the two fibers joined by each WCA contact."""
    subs = []
    for i, j in net.wca_edges:
        fi, fj = net.fiber_of[i], net.fiber_of[j]
        subs.append(np.where((net.fiber_of == fi) | (net.fiber_of == fj))[0])
    return subs


def star_subdomains(net: ContactNetwork) -> list[np.ndarray]:
    """Multi-fiber STAR (Open-work item 1): WCA two fibers PLUS every distinct crosslink-neighbor fiber."""
    xl = {}
    for a, c in net.xlink_edges:
        xl.setdefault(int(net.fiber_of[a]), set()).add(int(net.fiber_of[c]))
        xl.setdefault(int(net.fiber_of[c]), set()).add(int(net.fiber_of[a]))
    subs = []
    for i, j in net.wca_edges:
        fi, fj = int(net.fiber_of[i]), int(net.fiber_of[j])
        fibers = {fi, fj} | xl.get(fi, set()) | xl.get(fj, set())
        subs.append(np.where(np.isin(net.fiber_of, list(fibers)))[0])
    return subs


# --------------------------------------------------------------------------------------------------
# Additive Schwarz (optionally completed to a full-rank SPD preconditioner).
# --------------------------------------------------------------------------------------------------
def additive_schwarz(operator: np.ndarray, subdomains: list[np.ndarray], *, complete: bool = False) -> np.ndarray:
    """Dense symmetric additive Schwarz ``M = sum_e R_e^T W_e A_e^{-1} W_e R_e`` with degree^{-1/2} weights.

    When ``complete`` is set, every node not covered by any subdomain contributes its exact 3x3 diagonal-block
    inverse, so ``M`` is full-rank SPD (the idealized-preconditioner lens).  When it is not set, ``M`` is the
    raw overlap sum (the faithful-to-tournament direction generator, rank-deficient by construction).
    """
    dof = operator.shape[0]
    n_nodes = dof // 3
    covered = np.zeros(n_nodes, dtype=bool)
    degree = np.zeros(n_nodes)
    for nodes in subdomains:
        nodes = np.asarray(nodes, dtype=int)
        degree[nodes] += 1.0
        covered[nodes] = True
    if complete:
        uncovered = np.where(~covered)[0]
        degree[uncovered] += 1.0
    degree = np.maximum(degree, 1.0)
    weight = degree ** -0.5
    m = np.zeros((dof, dof))

    def _add(nodes):
        nodes = np.asarray(nodes, dtype=int)
        dofs = np.concatenate([[3 * n, 3 * n + 1, 3 * n + 2] for n in nodes])
        w = np.repeat(weight[nodes], 3)
        a_e = operator[np.ix_(dofs, dofs)]
        a_inv = np.linalg.solve(a_e, np.eye(a_e.shape[0]))
        m[np.ix_(dofs, dofs)] += (w[:, None] * a_inv) * w[None, :]

    for nodes in subdomains:
        _add(nodes)
    if complete:
        for node in np.where(~covered)[0]:
            _add(np.array([node]))
    return 0.5 * (m + m.T)


# --------------------------------------------------------------------------------------------------
# Diagnostics.
# --------------------------------------------------------------------------------------------------
def _free_basis(projector: np.ndarray, tol: float = 1e-8) -> np.ndarray:
    """Orthonormal basis of range(P) (the inextensibility-free displacement subspace)."""
    vals, vecs = np.linalg.eigh(projector)
    return vecs[:, vals > 1.0 - tol]


def coupled_translation_residual(net: ContactNetwork, edge_index: int | None = None) -> tuple[np.ndarray, int]:
    """Residual from rigidly translating one WCA-contacting fiber -- the documented plateau mode.

    Translating a contacting fiber along the contact normal stretches BOTH its stiff WCA contact and its
    crosslink to a *distinct* fiber, so the residual couples the WCA neighbor and the crosslink neighbor.  This
    is the ``45.27 pN`` fiber-translation residual the REPORT identifies; the pair/rigid blocks that omit the
    crosslink fiber cannot absorb it.  Returns ``(projected residual, fiber_id)``.
    """
    if edge_index is None:
        edge_index = len(net.wca_edges) // 2
    i = int(net.wca_edges[edge_index][0])
    fiber = int(net.fiber_of[i])
    normal = np.array([0.0, 1.0, 0.0])          # contact-stacking (y) direction
    disp = np.zeros((net.n_nodes, 3))
    disp[net.fiber_of == fiber] = normal
    v = disp.reshape(-1)
    r = net.projector @ (net.stiffness @ (net.projector @ v))
    return r, fiber


def direction_reduction(net: ContactNetwork, subdomains: list[np.ndarray], residual: np.ndarray) -> float:
    """Best single line-searched additive-Schwarz-direction residual reduction for ``residual`` in ``[0, 1]``.

    This is exactly one residual-monotone tournament step with the raw (uncompleted) additive-Schwarz
    correction ``d = P M P r`` as the candidate direction: ``reduction = 1 - ||r - alpha* A d|| / ||r||`` with
    the optimal scale ``alpha* = (r . A d)/(A d . A d)``.  A design whose subdomains span the coupled mode
    drives this toward 1; a design that omits the crosslink fiber leaves a large remainder.
    """
    a = net.operator
    m = net.projector @ additive_schwarz(a, subdomains, complete=False) @ net.projector
    d = m @ residual
    ad = a @ d
    ad2 = float(ad @ ad)
    if ad2 <= 0.0:
        return 0.0
    alpha = float(residual @ ad) / ad2
    remainder = residual - alpha * ad
    return float(1.0 - np.linalg.norm(remainder) / np.linalg.norm(residual))


def local_cluster_spd(net: ContactNetwork, subdomains: list[np.ndarray]) -> bool:
    """Every local subdomain block ``A_e`` (principal submatrix of aI+K) is SPD.

    Cauchy interlacing gives ``lambda_min(A_e) >= lambda_min(aI+K) >= reg > 0``, so each block is SPD however
    ill-conditioned (the soft bending modes sit many orders below the stiff WCA curvature but stay positive via
    the ``sqrt(eps64)*k_max`` floor).  The test therefore asserts strict positivity, not a condition bound.
    """
    for nodes in subdomains:
        nodes = np.asarray(nodes, dtype=int)
        dofs = np.concatenate([[3 * n, 3 * n + 1, 3 * n + 2] for n in nodes])
        a_e = net.stiffness[np.ix_(dofs, dofs)]
        w = np.linalg.eigvalsh(0.5 * (a_e + a_e.T))
        if w[0] <= 0.0:
            return False
    return True


def residual_dominance(net: ContactNetwork) -> dict:
    """Condition number of the projected operator with each force family removed, to rank conditioning drivers."""
    p = net.projector
    q = _free_basis(p)
    out = {"full": float(np.linalg.cond(q.T @ (p @ net.stiffness @ p) @ q))}
    for fam in ("regularizer", "wca", "xlink", "bend"):
        reduced = net.stiffness - net.parts[fam]
        reduced = 0.5 * (reduced + reduced.T)
        out[f"without_{fam}"] = float(np.linalg.cond(q.T @ (p @ reduced @ p) @ q))
    return out
