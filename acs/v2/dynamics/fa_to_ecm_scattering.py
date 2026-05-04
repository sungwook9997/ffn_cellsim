"""FA→ECM bilinear single-point scatter (Hard Blocker #3).

Implements the locked Hard Blocker #3 design from
``docs/v2_hard_blocker_3_fa_to_ecm_scattering_locked.md`` and the
pre-execution Sanity Gate
``docs/v2_fa_to_ecm_scattering_sanity_gate.md``: a pure function
that bilinearly scatters per-FA cell-on-substrate traction onto the
ECM grid as a vector traction-density field in ``nN/μm²``.

The function deliberately does **not**:

- mutate the input ECM state (pure function);
- call ``accumulate_prescribed_traction`` or any other accumulator
  / time-integration primitive (single-step deposit only);
- update any ECM response field
  (``stiffness_kpa`` / ``fiber_density`` / ``orientation_tensor``
  / ``accumulated_traction_nNs_per_um2`` / ``ligand_density``);
- scalarize the output (returns vector ``(nx, ny, 2)``, not a
  magnitude / norm field — the response law in Phase E owns the
  reduction choice per Hard Rule 11);
- silently clamp an out-of-grid FA to the nearest cell (raises
  ``fa_position_outside_ecm_grid`` instead);
- apply a Gaussian / disc / Boussinesq footprint (single-point
  bilinear only; footprint convolution deferred to Phase E).

Sanity Gate scope (acs/v2/dynamics/fa_to_ecm_scattering.py):

- §1 dimensional: per-FA per-cell ``[nN] · [dimensionless] /
  [μm²] = [nN/μm²]``. Single Rule 10 chain with no per-volume vs
  per-area mismatch. Spelled out in
  ``docs/v2_fa_to_ecm_scattering_sanity_gate.md`` §1.
- §2 boundary: every failure_kind from the Sanity Gate §2 is
  raised before any state computation; ``ecm.validate()`` runs at
  function entry per the local dynamics precedent
  (``acs.v2.dynamics.ecm_open_loop``).
- §3 conservation: component-wise vector conservation
  ``Σ_grid density · cell_area == Σ_FA traction`` to float64
  round-off. Pure-function structural tests on the input ECM and
  on a mocked ``accumulate_prescribed_traction`` enforce
  forbidden-side-effect constraints.
- §4 numerical: float64 throughout; per-FA work ``O(8)`` for
  interior bilinear; ``_FA_POSITION_BOUNDARY_TOL_UM`` is the only
  new numeric and is recomputed per call so it stays
  grid-invariant in relative terms.
- §5 sign: bilinear weights non-negative by construction;
  renormalized weights inherit; output preserves input traction
  sign.
- §6 measurement-protocol: vector output ``(nx, ny, 2)`` is the
  locked modality; no scalarization at scatter time. A runtime
  meta-test
  ``test_scatter_does_not_satisfy_closed_loop_response_law``
  enforces the Hard Rule 11 boundary at the test layer.

Magic-Number Block: only ``_FA_POSITION_BOUNDARY_TOL_UM = 1e-12 *
max(nx*dx, ny*dy)`` is new; documented as a relative float64
numerical tie-break (derivable, grid-invariant in relative terms,
not chosen to fit any test target).
"""

from __future__ import annotations

import math
from typing import Literal

import numpy as np

from acs.v2.ecm_substrate import ECMSubstrateState
from acs.v2.focal_adhesion import FocalAdhesionState

_BOUNDARY_TOL_RELATIVE = 1e-12

FAToECMScatteringFailureKind = Literal[
    "fa_position_outside_ecm_grid",
    "non_finite_fa_position",
    "non_finite_fa_traction",
]


class FAToECMScatteringError(ValueError):
    """Validation error raised by
    :func:`scatter_fa_traction_to_ecm_bilinear` when an FA's
    position or traction violates the locked scatter contract.
    Carries a machine-readable :attr:`failure_kind`."""

    def __init__(self, failure_kind: str, message: str) -> None:
        self.failure_kind = failure_kind
        super().__init__(message)


def _boundary_tol_um(ecm: ECMSubstrateState) -> float:
    """Recompute the numerical tie-break tolerance per ECM call so
    the inclusion check stays grid-invariant in relative terms.
    Documented in ``docs/v2_fa_to_ecm_scattering_sanity_gate.md`` §1."""

    nx, ny = ecm.grid_shape
    dx = float(ecm.spacing_um)
    return _BOUNDARY_TOL_RELATIVE * max(nx * dx, ny * dx)


def scatter_fa_traction_to_ecm_bilinear(
    adhesions: tuple[FocalAdhesionState, ...] | list[FocalAdhesionState],
    ecm: ECMSubstrateState,
) -> np.ndarray:
    """Bilinear single-point scatter of per-FA cell-on-substrate
    traction onto the ECM grid as a vector traction-density field.

    Returns an ``(nx, ny, 2)`` float64 array in ``nN/μm²``. The
    function is pure: ``ecm`` is not mutated, no accumulator is
    called, no response field is updated, and no scalarization is
    performed at this layer.

    Coordinate convention (locked design §2): cell-centered grid
    with ``x_i = origin_x + (i + 0.5) * dx``, inclusive physical
    footprint ``[origin_x, origin_x + nx*dx] × [origin_y,
    origin_y + ny*dy]``. Boundary half-cell positions use a
    truncated stencil with renormalized weights (locked design §3),
    preserving component-wise vector conservation at the boundary.

    Out-of-grid policy (locked design §4): a position outside the
    inclusive footprint plus the per-call numerical tie-break
    tolerance ``_FA_POSITION_BOUNDARY_TOL_UM`` raises
    ``FAToECMScatteringError(failure_kind=
    "fa_position_outside_ecm_grid")`` — no silent clamp.

    Failure modes:
    - ``non_finite_fa_position`` — any FA's ``position_um_xy`` is
      NaN / ±inf.
    - ``non_finite_fa_traction`` — any FA's
      ``traction_force_nN_xy`` is NaN / ±inf.
    - ``fa_position_outside_ecm_grid`` — any FA's position lies
      strictly outside the inclusive footprint plus tolerance.
    - Inherited ``ECMSubstrateState.validate()`` failures (e.g.,
      non-finite ``spacing_um``) propagate from the entry-point
      ``ecm.validate()`` call.

    Args:
        adhesions: tuple or list of :class:`FocalAdhesionState`.
        ecm: :class:`ECMSubstrateState`. Validated at entry.

    Returns:
        ``np.ndarray`` of shape ``(nx, ny, 2)`` and dtype
        ``np.float64``, vector traction density in ``nN/μm²`` per
        ECM cell. Component 0 is the x-component, component 1 is
        the y-component.
    """

    # Defensive ECM validation at the public dynamics boundary,
    # mirroring ``acs/v2/dynamics/ecm_open_loop.py`` precedent.
    ecm.validate()

    nx, ny = ecm.grid_shape
    dx = float(ecm.spacing_um)
    cell_area_um2 = dx * dx
    origin_x, origin_y = ecm.origin_um_xy
    origin_x = float(origin_x)
    origin_y = float(origin_y)
    tol = _boundary_tol_um(ecm)
    x_min = origin_x - tol
    x_max = origin_x + nx * dx + tol
    y_min = origin_y - tol
    y_max = origin_y + ny * dx + tol

    out = np.zeros((nx, ny, 2), dtype=np.float64)
    if not adhesions:
        return out

    for fa in adhesions:
        x_fa, y_fa = fa.position_um_xy
        x_fa = float(x_fa)
        y_fa = float(y_fa)
        if not (math.isfinite(x_fa) and math.isfinite(y_fa)):
            raise FAToECMScatteringError(
                "non_finite_fa_position",
                f"FA {fa.adhesion_id!r} position_um_xy={fa.position_um_xy!r} "
                f"is not finite",
            )
        f_x, f_y = fa.traction_force_nN_xy
        f_x = float(f_x)
        f_y = float(f_y)
        if not (math.isfinite(f_x) and math.isfinite(f_y)):
            raise FAToECMScatteringError(
                "non_finite_fa_traction",
                f"FA {fa.adhesion_id!r} traction_force_nN_xy="
                f"{fa.traction_force_nN_xy!r} is not finite",
            )
        if not (x_min <= x_fa <= x_max) or not (y_min <= y_fa <= y_max):
            raise FAToECMScatteringError(
                "fa_position_outside_ecm_grid",
                f"FA {fa.adhesion_id!r} at ({x_fa!r}, {y_fa!r}) is outside "
                f"ECM physical footprint [{origin_x}, {origin_x + nx*dx}] × "
                f"[{origin_y}, {origin_y + ny*dx}] (tolerance {tol!r})",
            )

        # Bracketing center indices: center i has x_i = origin_x +
        # (i + 0.5) * dx, so the lower bracket index is
        # i_lo = floor((x_fa - origin_x) / dx - 0.5).
        u = (x_fa - origin_x) / dx - 0.5
        v = (y_fa - origin_y) / dx - 0.5
        i_lo = int(math.floor(u))
        j_lo = int(math.floor(v))
        fx = u - i_lo
        fy = v - j_lo

        # Raw bilinear stencil with weights summing to 1 in the
        # interior. Some (i, j) pairs may be outside [0, nx) ×
        # [0, ny) when the FA is in a boundary half-cell.
        candidates = (
            (i_lo, j_lo, (1.0 - fx) * (1.0 - fy)),
            (i_lo + 1, j_lo, fx * (1.0 - fy)),
            (i_lo, j_lo + 1, (1.0 - fx) * fy),
            (i_lo + 1, j_lo + 1, fx * fy),
        )
        valid: list[tuple[int, int, float]] = [
            (k, l, w)
            for (k, l, w) in candidates
            if 0 <= k < nx and 0 <= l < ny
        ]
        if not valid:
            # Should be unreachable inside the inclusive footprint
            # + tolerance. If hit, treat as out-of-grid (defensive).
            raise FAToECMScatteringError(
                "fa_position_outside_ecm_grid",
                f"FA {fa.adhesion_id!r} at ({x_fa!r}, {y_fa!r}) bilinear "
                f"stencil has no in-grid centers; this should be unreachable "
                f"within the inclusive footprint plus tolerance",
            )

        weight_sum = sum(w for (_k, _l, w) in valid)
        # The interior case has weight_sum == 1.0 by construction;
        # the boundary half-cell case has weight_sum < 1.0 and is
        # renormalized so the deposited weights always sum to 1.0
        # (preserves component-wise vector conservation).
        for k, l, w in valid:
            w_norm = w / weight_sum
            out[k, l, 0] += w_norm * f_x / cell_area_um2
            out[k, l, 1] += w_norm * f_y / cell_area_um2

    return out
