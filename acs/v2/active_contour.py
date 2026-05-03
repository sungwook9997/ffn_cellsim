"""Active contour state and parameters for v2 P1 alpha (cortex + area, no dynamics).

Schema and analytic stability bound only. No physics force computation lives
here; that goes in :mod:`acs.v2.dynamics.active_contour`. This split keeps
the parameter contract and the rate bound together (they share the lock /
gate dimensional analysis) while leaving the actual time-stepping to the
dynamics module.

Reference: ``docs/v2_p1_derivation_locked.md`` §1 §2 §3 §4 (lock) and
``docs/v2_p1_active_contour_sanity_gate.md`` §1 §2 §3 §6 §7 (gate).

Sanity Gate scope (active_contour.py):

- §1 dimensional analysis: ``lambda_c_resolved_nN`` in nN,
  ``xi_line_resolved_nN_s_per_um2`` in nN·s/μm², ``k_a_nN_per_um``
  in nN/μm. Reductions verified in :func:`compute_rate_max`:
  Jacobian rows in [nN/μm], drag in [nN·s/μm], rate in [1/s].
- §2 timestep stability: :func:`compute_rate_max` returns the analytic
  Gershgorin row-sum upper bound and reports ``dt_cell · rate_max``
  against the named ``_DT_RATE_SAFETY_MARGIN``.
- §3 boundary cases: parameters validate cover under/over-specification
  of the dual-path cortex/drag config, non-positive ``k_a``, target
  area, ``dt_cell``, and ``n_vertices < 3``.
- §6 measurement-protocol consistency: :meth:`ActiveContourState.to_measurement_boundary`
  routes through :class:`MeasurementBoundary.from_array` so the canonical
  world_um_y_up CCW contract is preserved.

Magic-Number Block:

- ``_DT_RATE_SAFETY_MARGIN = 0.5`` is a named pre-run nonlinear safety
  margin: factor of 2 below the energy-monotone sufficient condition
  ``dt · rate_max ≤ 1`` (which is itself stricter than the linearised
  stability bound ``dt · rate_max < 2``). Not chosen to make any test
  pass (gate §2, §7).
- ``_AGREEMENT_TOLERANCE_RELATIVE = 1e-9`` is a named numerical
  validation tie-break for the float64 round-off in derived-vs-direct
  parameter equality checks (lock §1, §3 dual-path provenance). Not a
  physics tunable.
- ``min_edge_um`` defaults to the same value as the
  :class:`MeasurementBoundary` schema floor; it is inherited, not a
  fresh physics default.

No other tunable numerics live in this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from acs.v2.measurement_boundary import MeasurementBoundary

_DT_RATE_SAFETY_MARGIN = 0.5
_AGREEMENT_TOLERANCE_RELATIVE = 1e-9


class ActiveContourParametersError(ValueError):
    """Validation error carrying a machine-readable failure kind."""

    def __init__(self, failure_kind: str, message: str) -> None:
        self.failure_kind = failure_kind
        super().__init__(message)


@dataclass(frozen=True, slots=True)
class ActiveContourParameters:
    """Locked v2 P1 alpha active-contour parameters per
    ``docs/v2_p1_derivation_locked.md`` §7.

    Cortex coefficient ``lambda_c`` is supplied either directly
    (``lambda_c_nN``) or as the derived product
    ``sigma_c_nN_per_um · effective_height_um``. Drag density
    ``xi_line`` is supplied either directly
    (``xi_line_nN_s_per_um2``) or as the derived product
    ``xi_areal_nN_s_per_um3 · effective_height_um``. If both paths
    are supplied for either coefficient, the direct and derived
    values must agree to within
    ``_AGREEMENT_TOLERANCE_RELATIVE`` of their average; otherwise
    validation fails with ``cortex_overspecified_disagreement`` or
    ``drag_overspecified_disagreement``.

    A single canonical ``effective_height_um`` is shared between
    cortex and drag derived paths so an ``ActiveContourParameters``
    instance never resolves to two different cell heights for the
    two coefficients.
    """

    k_a_nN_per_um: float
    target_area_um2: float
    dt_cell_s: float
    n_vertices: int
    lambda_c_nN: Optional[float] = None
    sigma_c_nN_per_um: Optional[float] = None
    xi_line_nN_s_per_um2: Optional[float] = None
    xi_areal_nN_s_per_um3: Optional[float] = None
    effective_height_um: Optional[float] = None
    min_edge_um: float = 1e-3

    def validate(self) -> "ActiveContourParameters":
        # Cortex and area coefficients are allowed to be exactly zero so the
        # locked Tests 1/2/3 (zero-force, area-only, cortex-only) can be
        # configured without a pseudo-zero magic value (lock §6, gate §4).
        # `target_area`, `dt_cell`, drag, `effective_height`, and `min_edge`
        # remain strictly positive: a zero target area is degenerate, a zero
        # timestep does no integration, and zero drag would make the
        # overdamped equation of motion undefined.
        if not _is_finite_non_negative(self.k_a_nN_per_um):
            raise ActiveContourParametersError(
                "k_a_underspecified",
                "k_a_nN_per_um must be a finite non-negative float",
            )
        if not _is_finite_positive(self.target_area_um2):
            raise ActiveContourParametersError(
                "target_area_underspecified",
                "target_area_um2 must be a finite positive float",
            )
        if not _is_finite_positive(self.dt_cell_s):
            raise ActiveContourParametersError(
                "dt_underspecified",
                "dt_cell_s must be a finite positive float",
            )
        if (
            isinstance(self.n_vertices, bool)
            or not isinstance(self.n_vertices, int)
            or self.n_vertices < 3
        ):
            raise ActiveContourParametersError(
                "n_vertices_underspecified",
                "n_vertices must be an integer ≥ 3",
            )
        if not _is_finite_positive(self.min_edge_um):
            raise ActiveContourParametersError(
                "min_edge_param",
                "min_edge_um must be a finite positive float",
            )

        cortex_direct = self.lambda_c_nN is not None
        cortex_derived = (
            self.sigma_c_nN_per_um is not None
            and self.effective_height_um is not None
        )
        if cortex_direct and not _is_finite_non_negative(self.lambda_c_nN):
            raise ActiveContourParametersError(
                "cortex_underspecified",
                "lambda_c_nN must be a finite non-negative float when supplied",
            )
        if self.sigma_c_nN_per_um is not None and not _is_finite_non_negative(
            self.sigma_c_nN_per_um
        ):
            raise ActiveContourParametersError(
                "cortex_underspecified",
                "sigma_c_nN_per_um must be a finite non-negative float when supplied",
            )
        if not (cortex_direct or cortex_derived):
            raise ActiveContourParametersError(
                "cortex_underspecified",
                "supply lambda_c_nN OR (sigma_c_nN_per_um + effective_height_um)",
            )
        if cortex_direct and cortex_derived:
            derived = self.sigma_c_nN_per_um * self.effective_height_um
            if not _agree_relative(self.lambda_c_nN, derived):
                raise ActiveContourParametersError(
                    "cortex_overspecified_disagreement",
                    f"lambda_c_nN={self.lambda_c_nN!r} disagrees with "
                    f"sigma_c·H_eff={derived!r}",
                )

        drag_direct = self.xi_line_nN_s_per_um2 is not None
        drag_derived = (
            self.xi_areal_nN_s_per_um3 is not None
            and self.effective_height_um is not None
        )
        if drag_direct and not _is_finite_positive(self.xi_line_nN_s_per_um2):
            raise ActiveContourParametersError(
                "drag_underspecified",
                "xi_line_nN_s_per_um2 must be a finite positive float when supplied",
            )
        if self.xi_areal_nN_s_per_um3 is not None and not _is_finite_positive(
            self.xi_areal_nN_s_per_um3
        ):
            raise ActiveContourParametersError(
                "drag_underspecified",
                "xi_areal_nN_s_per_um3 must be a finite positive float when supplied",
            )
        if not (drag_direct or drag_derived):
            raise ActiveContourParametersError(
                "drag_underspecified",
                "supply xi_line_nN_s_per_um2 OR "
                "(xi_areal_nN_s_per_um3 + effective_height_um)",
            )
        if drag_direct and drag_derived:
            derived = self.xi_areal_nN_s_per_um3 * self.effective_height_um
            if not _agree_relative(self.xi_line_nN_s_per_um2, derived):
                raise ActiveContourParametersError(
                    "drag_overspecified_disagreement",
                    f"xi_line_nN_s_per_um2={self.xi_line_nN_s_per_um2!r} "
                    f"disagrees with xi_areal·H_eff={derived!r}",
                )

        if (cortex_derived or drag_derived) and not _is_finite_positive(
            self.effective_height_um
        ):
            raise ActiveContourParametersError(
                "effective_height_underspecified",
                "effective_height_um must be a finite positive float when "
                "any derived path is used",
            )
        return self

    @property
    def lambda_c_resolved_nN(self) -> float:
        if self.lambda_c_nN is not None:
            return float(self.lambda_c_nN)
        return float(self.sigma_c_nN_per_um) * float(self.effective_height_um)

    @property
    def xi_line_resolved_nN_s_per_um2(self) -> float:
        if self.xi_line_nN_s_per_um2 is not None:
            return float(self.xi_line_nN_s_per_um2)
        return float(self.xi_areal_nN_s_per_um3) * float(self.effective_height_um)


@dataclass
class ActiveContourState:
    """Mutable simulation state for the active contour.

    Vertices are stored as a writable ``(N, 2)`` float64 array in
    canonical world_um_y_up coordinates. The state never exposes its
    raw vertices to downstream measurement code; downstream calls
    :meth:`to_measurement_boundary` and consumes the immutable
    :class:`MeasurementBoundary`.

    ``cell_id`` is carried so artifact writers can tag the polygon
    when the active contour is embedded in a single-cell or cluster
    context.
    """

    cell_id: str
    vertices_xy_um: np.ndarray
    params: ActiveContourParameters
    step_count: int = 0

    def __post_init__(self) -> None:
        # Reject embedding an unvalidated parameters object: state must
        # never silently carry an invalid contract into compute_rate_max
        # or downstream dynamics.
        self.params.validate()
        if not isinstance(self.cell_id, str) or not self.cell_id.strip():
            raise ActiveContourParametersError(
                "cell_id_underspecified",
                "cell_id must be a non-empty string",
            )
        if (
            isinstance(self.step_count, bool)
            or not isinstance(self.step_count, int)
            or self.step_count < 0
        ):
            raise ActiveContourParametersError(
                "step_count_underspecified",
                "step_count must be a non-negative integer",
            )
        verts = np.asarray(self.vertices_xy_um, dtype=np.float64)
        if verts.ndim != 2 or verts.shape[1] != 2:
            raise ActiveContourParametersError(
                "vertices_shape",
                f"vertices_xy_um must have shape (N, 2), got {verts.shape!r}",
            )
        if verts.shape[0] != self.params.n_vertices:
            raise ActiveContourParametersError(
                "vertices_count_mismatch",
                f"vertices_xy_um has {verts.shape[0]} rows, "
                f"params.n_vertices={self.params.n_vertices}",
            )
        if not np.isfinite(verts).all():
            raise ActiveContourParametersError(
                "non_finite_vertices",
                "vertices_xy_um must contain only finite values",
            )
        # Reseat as a writable array so callers cannot smuggle a
        # frozen view.
        self.vertices_xy_um = np.array(verts, dtype=np.float64, copy=True)
        # Pre/post-step canonical contract: every ActiveContourState is
        # always a valid MeasurementBoundary polygon (Hard Rule 11).
        # A state instantiated from a self-intersecting / CW / pinched
        # input therefore raises MeasurementBoundaryError immediately at
        # construction time, before any force evaluation can run.
        self.to_measurement_boundary()

    def to_measurement_boundary(
        self, *, source_modality: str = "active_contour"
    ) -> MeasurementBoundary:
        return MeasurementBoundary.from_array(
            self.vertices_xy_um,
            coordinate_convention="world_um_y_up",
            source_modality=source_modality,
            min_edge_um=self.params.min_edge_um,
        )

    def area_um2(self) -> float:
        x = self.vertices_xy_um[:, 0]
        y = self.vertices_xy_um[:, 1]
        return float(0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))

    def perimeter_um(self) -> float:
        diffs = np.roll(self.vertices_xy_um, -1, axis=0) - self.vertices_xy_um
        return float(np.linalg.norm(diffs, axis=1).sum())

    def edge_lengths_um(self) -> np.ndarray:
        diffs = np.roll(self.vertices_xy_um, -1, axis=0) - self.vertices_xy_um
        return np.linalg.norm(diffs, axis=1)


def compute_rate_max(state: ActiveContourState) -> dict:
    """Return the analytic Gershgorin upper bound on ``rate_max`` and
    the runtime contract diagnostics required by the Sanity Gate.

    Per gate §2 and lock §4 the cortex and area force-Jacobian rows are
    computed in ``[nN/μm]`` and then divided by per-vertex drag
    ``ζ_i [nN·s/μm]`` exactly once to give the rate row in ``[1/s]``.
    The returned ``rate_max`` is the maximum rate row.

    The returned dict contains:

    - ``rate_max_s_inv``: scalar [1/s]
    - ``dt_cell_s``: scalar [s]
    - ``dt_rate_product``: ``dt · rate_max`` (dimensionless)
    - ``margin_ratio``: ``dt · rate_max / _DT_RATE_SAFETY_MARGIN``
      (must be ≤ 1 for the runtime contract; not enforced here).
    - ``jacobian_row_cortex_nN_per_um``: per-vertex array [nN/μm],
      summed BEFORE dividing by drag.
    - ``jacobian_row_area_nN_per_um``: per-vertex array [nN/μm],
      summed BEFORE dividing by drag.
    - ``jacobian_row_total_nN_per_um``: ``cortex + area`` row sum,
      again BEFORE the single division by drag.
    - ``zeta_nN_s_per_um``: per-vertex drag array [nN·s/μm].
    - ``rate_per_vertex_s_inv``: per-vertex rate array [1/s].
    - ``min_edge_length_um``: scalar [μm].
    """

    params = state.params
    verts = state.vertices_xy_um
    n = params.n_vertices

    # Edge lengths. ell_edge[i] = |r_{i+1} - r_i|.
    diffs = np.roll(verts, -1, axis=0) - verts
    ell_edge = np.linalg.norm(diffs, axis=1)
    if (ell_edge <= 0.0).any():
        raise ActiveContourParametersError(
            "edge_length_zero",
            "active contour state has zero-length edges; "
            "cannot bound rate_max",
        )
    # Vertex control length: half-sum of adjacent edges.
    # ell_v[i] = 0.5 * (ell_edge[i] + ell_edge[i-1]).
    ell_v = 0.5 * (ell_edge + np.roll(ell_edge, 1))

    # Cortex Jacobian row bound (lock §4): J_i^c <= 2 lambda_c (1/ell_i + 1/ell_{i-1}).
    lambda_c = params.lambda_c_resolved_nN
    inv_ell_i = 1.0 / ell_edge
    inv_ell_im1 = 1.0 / np.roll(ell_edge, 1)
    jacobian_row_cortex = 2.0 * lambda_c * (inv_ell_i + inv_ell_im1)

    # Area Jacobian row exact (lock §4):
    # J_i^A = (K_A/A_0) |g_i| sum_j |g_j| + |p_A| * H^A_row_i.
    # H^A_row_i = 1.0 (closed-form: blocks H^A_{i, i-1} and H^A_{i, i+1}
    # each have operator norm 0.5; H^A_{i,i} = 0; H^A_{i,j} = 0 for
    # |i-j| > 1).
    k_a = float(params.k_a_nN_per_um)
    target_area = float(params.target_area_um2)
    current_area = state.area_um2()
    p_a = k_a * (current_area - target_area) / target_area  # [nN/μm]
    grad = np.column_stack(
        (
            0.5 * (np.roll(verts[:, 1], -1) - np.roll(verts[:, 1], 1)),
            0.5 * (np.roll(verts[:, 0], 1) - np.roll(verts[:, 0], -1)),
        )
    )  # [μm]
    grad_norm = np.linalg.norm(grad, axis=1)
    sum_grad_norms = float(grad_norm.sum())
    jacobian_row_area = (k_a / target_area) * grad_norm * sum_grad_norms + abs(p_a) * 1.0

    jacobian_row_total = jacobian_row_cortex + jacobian_row_area
    zeta = params.xi_line_resolved_nN_s_per_um2 * ell_v
    rate_per_vertex = jacobian_row_total / zeta
    rate_max = float(rate_per_vertex.max())

    dt_cell = float(params.dt_cell_s)
    dt_rate_product = dt_cell * rate_max
    return {
        "rate_max_s_inv": rate_max,
        "dt_cell_s": dt_cell,
        "dt_rate_product": dt_rate_product,
        "margin_ratio": dt_rate_product / _DT_RATE_SAFETY_MARGIN,
        "jacobian_row_cortex_nN_per_um": jacobian_row_cortex,
        "jacobian_row_area_nN_per_um": jacobian_row_area,
        "jacobian_row_total_nN_per_um": jacobian_row_total,
        "zeta_nN_s_per_um": zeta,
        "rate_per_vertex_s_inv": rate_per_vertex,
        "min_edge_length_um": float(ell_edge.min()),
    }


def regular_polygon_vertices(n: int, radius_um: float, *, center=(0.0, 0.0)) -> np.ndarray:
    """Build a regular ``n``-gon at ``radius_um`` for tests / harness.

    Returns CCW vertices in canonical world_um_y_up. Phase is chosen so
    the polygon is not aligned with axes (avoids accidental symmetry
    cancellations in the area gradient test).
    """

    if n < 3:
        raise ActiveContourParametersError(
            "n_vertices_underspecified",
            "regular polygon requires n >= 3",
        )
    if not _is_finite_positive(radius_um):
        raise ActiveContourParametersError(
            "radius_underspecified",
            "radius_um must be a finite positive float",
        )
    angles = 2.0 * np.pi * np.arange(n, dtype=np.float64) / n + np.pi / (2 * n)
    cx, cy = float(center[0]), float(center[1])
    return np.column_stack(
        (cx + radius_um * np.cos(angles), cy + radius_um * np.sin(angles))
    )


def ellipse_polygon_vertices(
    n: int, semi_a_um: float, semi_b_um: float, *, center=(0.0, 0.0)
) -> np.ndarray:
    """Build an ellipse-sampled polygon for the Test 4 diagnostic."""

    if n < 3:
        raise ActiveContourParametersError(
            "n_vertices_underspecified",
            "ellipse polygon requires n >= 3",
        )
    if not (_is_finite_positive(semi_a_um) and _is_finite_positive(semi_b_um)):
        raise ActiveContourParametersError(
            "radius_underspecified",
            "semi_a_um and semi_b_um must be finite positive",
        )
    angles = 2.0 * np.pi * np.arange(n, dtype=np.float64) / n
    cx, cy = float(center[0]), float(center[1])
    return np.column_stack(
        (cx + semi_a_um * np.cos(angles), cy + semi_b_um * np.sin(angles))
    )


def _is_finite_positive(value) -> bool:
    if value is None:
        return False
    if isinstance(value, bool):
        return False
    try:
        f = float(value)
    except (TypeError, ValueError):
        return False
    return np.isfinite(f) and f > 0.0


def _is_finite_non_negative(value) -> bool:
    if value is None:
        return False
    if isinstance(value, bool):
        return False
    try:
        f = float(value)
    except (TypeError, ValueError):
        return False
    return np.isfinite(f) and f >= 0.0


def _agree_relative(a: float, b: float) -> bool:
    a = float(a)
    b = float(b)
    if not (np.isfinite(a) and np.isfinite(b)):
        return False
    avg = 0.5 * (abs(a) + abs(b))
    if avg == 0.0:
        return abs(a - b) <= _AGREEMENT_TOLERANCE_RELATIVE
    return abs(a - b) <= _AGREEMENT_TOLERANCE_RELATIVE * avg
