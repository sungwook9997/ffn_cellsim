"""ECM substrate state schema.

Canonical Phase 1 ``ECMSubstrateState`` per
``docs/v2_phase1_plan_consolidated.md`` §5.3. Phase 1 represents the ECM
as a coarse continuous field on a regular grid. Full collagen-fiber
graph representations are deferred.

Dynamics (stiffening, alignment, density change in response to traction)
are explicitly out of P0 scope. The ECM Sanity Gate must derive the
unit chain between kPa, pN/um^2 (force per area), traction density, and
grid-cell storage before any executable update is added (Plan §6.4 and
Hard Rule 10). This schema only validates state identity and the
numerical invariants of the stored arrays.

Sanity Gate scope: non-physics schema module. Full physics 6-item gate
is N/A; this module only owns the boundary-case checks (positive grid
spacing, matching shapes, finiteness, density bounds, orientation
tensor symmetry, non-negative scalar traction history) and the
measurement-protocol consistency that the orientation tensor is stored
per grid cell (shape ``(nx, ny, 2, 2)``) in dimensionless units.

Magic-Number Block: this module declares no tunable numeric. N/A.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

# Numerical validation tolerances. Both are numerical tie-breaks for
# float64 round-off in the schema layer, NOT physical or fitted tunables.
_ORIENTATION_SYMMETRY_TOL = 1e-12  # max |T - T.T| allowed at any grid cell
_ORIENTATION_BOUND = 1.0  # |T_ij| <= 1.0 for the dimensionless orientation tensor

ECMSource = Literal["synthetic", "image_derived", "literature_default"]


@dataclass(frozen=True)
class ECMSubstrateState:
    """Schema for the coarse continuous ECM field on a regular xy grid."""

    origin_um_xy: tuple[float, float]
    spacing_um: float
    stiffness_kpa: NDArray[np.float64]
    ligand_density: NDArray[np.float64]
    fiber_density: NDArray[np.float64]
    orientation_tensor: NDArray[np.float64]
    accumulated_traction_nNs_per_um2: NDArray[np.float64]
    source: ECMSource = "synthetic"

    def validate(self) -> None:
        if len(self.origin_um_xy) != 2 or not all(
            np.isfinite(v) for v in self.origin_um_xy
        ):
            raise ValueError("origin_um_xy must contain exactly 2 finite values")
        if not np.isfinite(self.spacing_um) or self.spacing_um <= 0.0:
            raise ValueError("spacing_um must be positive and finite")

        scalar_arrays = {
            "stiffness_kpa": self.stiffness_kpa,
            "ligand_density": self.ligand_density,
            "fiber_density": self.fiber_density,
            "accumulated_traction_nNs_per_um2": self.accumulated_traction_nNs_per_um2,
        }
        shapes = {name: np.asarray(arr).shape for name, arr in scalar_arrays.items()}
        unique_shapes = set(shapes.values())
        if len(unique_shapes) != 1:
            raise ValueError(f"ECM scalar grids must share shape, got {shapes!r}")
        grid_shape = unique_shapes.pop()
        if len(grid_shape) != 2 or grid_shape[0] <= 0 or grid_shape[1] <= 0:
            raise ValueError(
                f"ECM scalar grids must have shape (nx>0, ny>0), got {grid_shape!r}"
            )

        for name, arr in scalar_arrays.items():
            arr_np = np.asarray(arr)
            if not np.isfinite(arr_np).all():
                raise ValueError(f"{name} must be finite at every grid cell")

        if (np.asarray(self.stiffness_kpa) < 0.0).any():
            raise ValueError("stiffness_kpa must be non-negative")
        ligand = np.asarray(self.ligand_density)
        if (ligand < 0.0).any() or (ligand > 1.0).any():
            raise ValueError("ligand_density must be in [0, 1]")
        fiber = np.asarray(self.fiber_density)
        if (fiber < 0.0).any() or (fiber > 1.0).any():
            raise ValueError("fiber_density must be in [0, 1]")
        if (np.asarray(self.accumulated_traction_nNs_per_um2) < 0.0).any():
            raise ValueError(
                "accumulated_traction_nNs_per_um2 must be non-negative when stored "
                "as a scalar history"
            )

        orientation = np.asarray(self.orientation_tensor)
        if orientation.shape != (*grid_shape, 2, 2):
            raise ValueError(
                f"orientation_tensor must have shape {(*grid_shape, 2, 2)!r}, "
                f"got {orientation.shape!r}"
            )
        if not np.isfinite(orientation).all():
            raise ValueError("orientation_tensor must be finite at every grid cell")
        symmetry_residual = np.max(
            np.abs(orientation - np.swapaxes(orientation, -1, -2))
        )
        if symmetry_residual > _ORIENTATION_SYMMETRY_TOL:
            raise ValueError(
                f"orientation_tensor must be symmetric (max |T - T.T| = {symmetry_residual!r})"
            )
        max_abs = float(np.max(np.abs(orientation)))
        if max_abs > _ORIENTATION_BOUND:
            raise ValueError(
                f"orientation_tensor components must satisfy |T_ij| <= "
                f"{_ORIENTATION_BOUND} (got max {max_abs!r}); the orientation "
                f"tensor is dimensionless and bounded"
            )

        if self.source not in {"synthetic", "image_derived", "literature_default"}:
            raise ValueError(
                "source must be synthetic, image_derived, or literature_default"
            )

    @property
    def grid_shape(self) -> tuple[int, int]:
        shape = np.asarray(self.stiffness_kpa).shape
        return (int(shape[0]), int(shape[1]))
