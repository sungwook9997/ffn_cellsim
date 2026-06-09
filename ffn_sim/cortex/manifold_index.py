"""Manifold broad-phase INDEX over the cortex surface (Option-A, search/coordinate).

H.7 2026-06-09 (PI: option A — wire the Layer-2 geometry manifold into the cortex as a
SEARCH/MEMORY/COORDINATE accelerator ONLY). This is the *index half*: a thin, persistent
wrapper around :class:`ffn_sim.cortex.surface_manifold.SurfaceManifold` that replaces the
per-tick global ``scipy.spatial.cKDTree`` rebuild every binder currently does with a
persistent bead→patch map + geodesic k-ring candidate gather, refreshed lazily.

SCOPE (ratified by the 2026-06-09 (a)/(b) adversarial study,
``H7_2D_MESH_AB_DECISION_2026-06-09.md``):

* This carries **NO mechanics** — no forces, no tension, no γ DOF, mesh edges are NOT
  filaments. It is the sanctioned Layer-2 use (``H7_SURFACE_MANIFOLD_EXPLICIT_CORTEX``).
* It is **NOT a γ fix** — the active-γ transmission floor is set by the soft scalar
  couplings (k_head_actin/k_intra), orthogonal to where/how candidates are found. This
  index changes the SEARCH and the shared coordinate frame, nothing in the force budget.
* The candidate set it returns is **byte-identical** to the global cKDTree within-reach
  query AT REFRESH TIME (broad-phase k-ring covers the reach ball, narrow-phase filters by
  the same Euclidean reach) — proven by ``h7_manifold_index_validate`` (VG-6).
* ⚠️ The CORRECTNESS HOLE the audit surfaced: a STALE ``bead_tri`` (lazy refresh) silently
  drops within-reach pairs once beads drift ~1 patch. This module therefore exposes
  :meth:`max_drift_since_refresh` and a :meth:`needs_refresh` policy, and the validation
  derives the safe refresh cadence. A binder MUST refresh before the drift gate trips.
* Performance honesty: at the ×40 mesoscale on CPU the Python k-ring gather is SLOWER than
  scipy's C cKDTree (measured); the win is conditional on GPU-residency. This module keeps
  the gather numpy-vectorised and a cupy port is the follow-on. Adopt for spatiality /
  shared coordinate frame first; speed only once GPU-resident.

Sanity Gate (CLAUDE.md): dimensional (lengths [m] in, dimensionless indices out; no force);
boundary (empty cloud → empty; k=0 → home patch only); conservation (none — geometry only);
measurement consistency (within-reach pair set == global query at refresh time, VG-6);
no-magic-number (k is DERIVED via ``kring_for_reach``; refresh threshold is a fraction of the
mesh circumradius, geometry-derived, not tuned).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ffn_sim.cortex.surface_manifold import SurfaceManifold


@dataclass
class ManifoldIndex:
    """Persistent bead→patch broad-phase index over a :class:`SurfaceManifold`.

    Attributes:
        manifold: the surface manifold (geometry only).
        reach: the physical candidate-search radius [m] (e.g. ``max_bind_dist`` or the
            ``√(A/n)`` bridge reach). Sets the derived k-ring hop count.
        safety_rings: extra k-ring rings for ball-coverage robustness (passed through to
            :meth:`SurfaceManifold.kring_for_reach`).
        refresh_drift_frac: refresh when the max bead drift since the last refresh exceeds
            ``refresh_drift_frac · max_circumradius`` (geometry-derived staleness gate; a
            bead that has drifted a full patch-circumradius may have changed home patch).
    """

    manifold: SurfaceManifold
    reach: float
    safety_rings: int = 1
    refresh_drift_frac: float = 0.5

    k: int = field(init=False)
    _bead_tri: np.ndarray | None = field(default=None, init=False, repr=False)
    _ref_pos: np.ndarray | None = field(default=None, init=False, repr=False)
    # per-patch k-ring (list of triangle arrays), cached once per manifold geometry:
    _patch_kring: list[np.ndarray] | None = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        if not (np.isfinite(self.reach) and self.reach > 0.0):
            raise ValueError(f"reach must be finite > 0; got {self.reach!r}")
        self.k = self.manifold.kring_for_reach(self.reach, self.safety_rings)

    # -- refresh / staleness ------------------------------------------------
    def refresh(self, pos: np.ndarray, *, fit: bool = False) -> "ManifoldIndex":
        """Rebuild the bead→patch map from the current bead cloud.

        Args:
            pos: bead positions, shape ``(n, 3)`` [m].
            fit: if True, first slave the manifold geometry to the cloud
                (:meth:`SurfaceManifold.fit_to_cloud`) — the deformable-shell update for a
                non-spherical / spreading cell. Default False (rigid sphere shell).
        """
        pos = np.asarray(pos, dtype=np.float64)
        if pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError(f"pos must be (n, 3); got {pos.shape}")
        if fit and pos.shape[0] > 0:
            self.manifold.fit_to_cloud(pos)
            self._patch_kring = None  # geometry moved → adjacency same but recache lazily
        self._bead_tri = self.manifold.nearest_patch(pos)
        self._ref_pos = pos.copy()
        return self

    def max_drift_since_refresh(self, pos: np.ndarray) -> float:
        """Max bead displacement [m] since the last :meth:`refresh` (the staleness signal)."""
        if self._ref_pos is None:
            return float("inf")
        pos = np.asarray(pos, dtype=np.float64)
        if pos.shape != self._ref_pos.shape:
            return float("inf")
        if pos.shape[0] == 0:
            return 0.0
        return float(np.max(np.linalg.norm(pos - self._ref_pos, axis=1)))

    @property
    def drift_threshold(self) -> float:
        """Geometry-derived staleness threshold [m] = frac · patch circumradius."""
        return self.refresh_drift_frac * float(self.manifold.max_circumradius)

    def needs_refresh(self, pos: np.ndarray) -> bool:
        """True if the bead cloud has drifted past the staleness gate (refresh required)."""
        return self.max_drift_since_refresh(pos) >= self.drift_threshold

    # -- candidate gather (broad phase) -------------------------------------
    def _ensure_patch_kring(self) -> None:
        if self._patch_kring is None:
            self._patch_kring = [
                self.manifold.patch_kring(t, self.k) for t in range(self.manifold.n_tri)
            ]

    def neighbor_pools(self) -> list[np.ndarray]:
        """Per-bead broad-phase candidate pool (bead indices homed in the home patch's
        k-ring), using the CURRENT (possibly stale) bead→patch map. Narrow-phase Euclidean
        filtering is the caller's job (this is broad-phase only)."""
        if self._bead_tri is None:
            raise RuntimeError("call refresh(pos) before neighbor_pools().")
        self._ensure_patch_kring()
        n = self._bead_tri.shape[0]
        # patch -> beads homed there
        patch_beads: list[list[int]] = [[] for _ in range(self.manifold.n_tri)]
        for b in range(n):
            patch_beads[int(self._bead_tri[b])].append(b)
        patch_beads_arr = [np.array(x, dtype=np.int64) for x in patch_beads]
        pools: list[np.ndarray] = []
        for b in range(n):
            tris = self._patch_kring[int(self._bead_tri[b])]
            parts = [patch_beads_arr[int(t)] for t in tris if patch_beads_arr[int(t)].size]
            pools.append(np.concatenate(parts) if parts else np.empty(0, dtype=np.int64))
        return pools

    def pairs_within_reach(
        self, pos: np.ndarray, *, reach: float | None = None, exclude_self: bool = True
    ) -> set[tuple[int, int]]:
        """Within-reach undirected pairs via the broad-phase k-ring pool + Euclidean
        narrow phase on ``pos``. This is the binder-facing query; at refresh time its
        result is identical to the global cKDTree within-reach query (VG-6).

        Uses the CURRENT bead→patch map — if stale, pairs can be MISSED (the audit hole);
        guard with :meth:`needs_refresh`.
        """
        r = float(self.reach if reach is None else reach)
        pos = np.asarray(pos, dtype=np.float64)
        pools = self.neighbor_pools()
        out: set[tuple[int, int]] = set()
        for i, pool in enumerate(pools):
            if pool.size == 0:
                continue
            d = np.linalg.norm(pos[pool] - pos[i], axis=1)
            for j in pool[d <= r]:
                j = int(j)
                if exclude_self and j == i:
                    continue
                out.add((i, j) if i < j else (j, i))
        return out


def press_onto_substrate(manifold: SurfaceManifold, z_basal: float) -> SurfaceManifold:
    """Deform a (suspended) sphere manifold into the ADHERENT cell shape: clamp every
    vertex below ``z_basal`` up onto the substrate plane z = z_basal, leaving the upper
    body rounded — a FLAT VENTRAL surface + rounded apical (PI 2026-06-09 adherent pivot).

    Uses only :meth:`SurfaceManifold.set_verts` (no edit to the shared manifold class) so
    the ventral region becomes a flat triangulated patch CONSISTENT with the flat ventral
    filament placement (removes the S1-curved / B1-flat mismatch that produced the spurious
    "curvature option"). GEOMETRY ONLY — no mechanics, no γ.

    Args:
        manifold: an icosphere (or any) surface manifold, modified in place.
        z_basal: the substrate plane height [m] (verts below it are clamped onto it).

    Returns:
        ``manifold`` (pressed in place via ``set_verts``).
    """
    verts = manifold.verts.copy()
    below = verts[:, 2] < z_basal
    verts[below, 2] = z_basal
    manifold.set_verts(verts)
    return manifold


def global_pairs_within_reach(
    pos: np.ndarray, reach: float, *, exclude_self: bool = True
) -> set[tuple[int, int]]:
    """Ground-truth within-reach undirected pairs via a global ``cKDTree`` (the query the
    index must reproduce at refresh time). Geometry only."""
    from scipy.spatial import cKDTree

    pos = np.asarray(pos, dtype=np.float64)
    if pos.shape[0] == 0:
        return set()
    tree = cKDTree(pos)
    out: set[tuple[int, int]] = set()
    for i, js in enumerate(tree.query_ball_point(pos, reach)):
        for j in js:
            if exclude_self and j == i:
                continue
            out.add((i, j) if i < j else (j, i))
    return out
