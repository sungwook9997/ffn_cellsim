"""Membrane--ERM--cortex preload capacity arithmetic (host acceptance algebra).

The production ``t0`` is a pressure-bearing equilibrium, not a force-free mesh.  Before spending CUDA time
on a relaxation, this module evaluates the algebra obtained from membrane tension plus an explicit bound-linker
population at a supplied per-link force basis.  It is a physiological upper bound **only if the caller supplies
a source-grounded single-ERM force**.  A continuum membrane tube-extraction force is not interchangeable with
that molecular input; the current Active Cell caller labels that historical substitution diagnostic and the
NG-3 gate rejects it.

This is acceptance algebra only.  It does not add an inward force, choose an ERM density, or replace explicit
tethers with a surface stress.  The runtime still obtains all membrane--cortex load transfer from individual
ERM spring states.

Sanity Gate:
    * Dimensions: ``2*gamma/R``, ``rho*f_rupt``, and ``pressure`` are pN/um^2 == Pa.  Counts are
      dimensionless and densities are um^-2.
    * Boundary cases: zero pressure requires zero tethers; membrane tension alone may close the load; zero
      tethers cannot support a positive residual pressure; non-positive geometry/rupture inputs are rejected.
    * Sign sense: pressure and required density increase together; membrane tension, rupture force, and
      explicit tether density increase capacity.
    * Grid invariance: the required *density* is independent of membrane triangulation area; only the integer
      count scales with the represented surface area.
    * No fitting: with a valid molecular input the result is a necessary capacity bound, never a production
      density prescription; without one it is diagnostic arithmetic only.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass

import numpy as np

__all__ = ["PreloadCapacity", "evaluate_preload_capacity", "triangle_surface_area"]


@dataclass(frozen=True)
class PreloadCapacity:
    """Load-capacity arithmetic whose production meaning depends on the supplied force provenance."""

    pressure_pa: float
    radius_um: float
    membrane_area_um2: float
    membrane_tension_pn_per_um: float
    membrane_laplace_pressure_pa: float
    residual_pressure_for_erm_pa: float
    erm_tether_count: int
    erm_density_per_um2: float
    erm_rupture_force_pn: float
    erm_pressure_capacity_upper_bound_pa: float
    total_pressure_capacity_upper_bound_pa: float
    pressure_capacity_ratio: float
    required_erm_density_min_per_um2: float
    required_erm_count_min: int
    mechanically_feasible_upper_bound: bool

    def ledger_fields(self) -> dict[str, float | int | bool]:
        """Flatten the immutable result for the production population ledger."""
        return {f"preload_{key}": value for key, value in asdict(self).items()}


def triangle_surface_area(verts: np.ndarray, faces: np.ndarray) -> float:
    """Return the exact flat-triangle area represented by a closed mesh [um^2]."""
    v = np.asarray(verts, dtype=np.float64)
    f = np.asarray(faces, dtype=np.int64)
    if v.ndim != 2 or v.shape[1] != 3:
        raise ValueError(f"verts must be (N, 3); got {v.shape}")
    if f.ndim != 2 or f.shape[1] != 3:
        raise ValueError(f"faces must be (M, 3); got {f.shape}")
    if f.size == 0:
        return 0.0
    area2 = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
    return float(0.5 * np.linalg.norm(area2, axis=1).sum())


def evaluate_preload_capacity(
    *,
    pressure_pa: float,
    radius_um: float,
    membrane_area_um2: float,
    membrane_tension_pn_per_um: float,
    erm_tether_count: int,
    erm_rupture_force_pn: float,
) -> PreloadCapacity:
    """Evaluate pressure-bearing capacity arithmetic for a supplied per-link force basis.

    The continuum membrane contributes ``2*gamma/R``.  The most optimistic possible explicit-ERM
    contribution is ``(N/A)*f_rupt``.  If and only if ``f_rupt`` is a source-grounded single-ERM value, passing
    this bound establishes necessary (not sufficient) feasibility and failing it proves impossibility.  The
    function cannot infer provenance from a scalar; callers must record and gate the force basis separately.
    """
    values = {
        "pressure_pa": pressure_pa,
        "radius_um": radius_um,
        "membrane_area_um2": membrane_area_um2,
        "membrane_tension_pn_per_um": membrane_tension_pn_per_um,
        "erm_rupture_force_pn": erm_rupture_force_pn,
    }
    if not all(np.isfinite(float(value)) for value in values.values()):
        raise ValueError(f"preload inputs must be finite; got {values}")
    if pressure_pa < 0.0:
        raise ValueError("pressure_pa must be nonnegative")
    if radius_um <= 0.0 or membrane_area_um2 <= 0.0 or erm_rupture_force_pn <= 0.0:
        raise ValueError("radius, membrane area, and ERM rupture force must be positive")
    if membrane_tension_pn_per_um < 0.0:
        raise ValueError("membrane tension must be nonnegative")
    if int(erm_tether_count) != erm_tether_count or erm_tether_count < 0:
        raise ValueError("erm_tether_count must be a nonnegative integer")

    count = int(erm_tether_count)
    p_mem = 2.0 * float(membrane_tension_pn_per_um) / float(radius_um)
    p_residual = max(float(pressure_pa) - p_mem, 0.0)
    density = count / float(membrane_area_um2)
    p_erm_upper = density * float(erm_rupture_force_pn)
    p_total_upper = p_mem + p_erm_upper
    required_density = p_residual / float(erm_rupture_force_pn)
    required_count = int(math.ceil(required_density * float(membrane_area_um2)))
    feasible = bool(count >= required_count)
    capacity_ratio = math.inf if pressure_pa == 0.0 else p_total_upper / float(pressure_pa)

    return PreloadCapacity(
        pressure_pa=float(pressure_pa),
        radius_um=float(radius_um),
        membrane_area_um2=float(membrane_area_um2),
        membrane_tension_pn_per_um=float(membrane_tension_pn_per_um),
        membrane_laplace_pressure_pa=p_mem,
        residual_pressure_for_erm_pa=p_residual,
        erm_tether_count=count,
        erm_density_per_um2=density,
        erm_rupture_force_pn=float(erm_rupture_force_pn),
        erm_pressure_capacity_upper_bound_pa=p_erm_upper,
        total_pressure_capacity_upper_bound_pa=p_total_upper,
        pressure_capacity_ratio=capacity_ratio,
        required_erm_density_min_per_um2=required_density,
        required_erm_count_min=required_count,
        mechanically_feasible_upper_bound=feasible,
    )
