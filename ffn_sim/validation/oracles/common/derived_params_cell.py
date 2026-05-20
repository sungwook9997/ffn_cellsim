"""Resolve KU-derived parameters for Phase 1 Unit 3 (single-cell cortex).

All formulas come from cited KUs; no fitting, no hand-tuning.

Geometry (KU-3.17):
    Cell circumference:        C = 2π R_cell                          [m]
    Fiber coverage ratio:      ρ_cov = N_f · L_f / C                  [dimensionless]
    Bead rest length:          ℓ₀ = L_f / (beads_per_fiber − 1)        [m]

Time integration (KU-3.16 quasi-static + KU-1.26 overdamped form):
    Per-bead drag γ_b is taken from ``cell.gamma_drag_per_bead`` (the
    Phase 1 lumped value; see the yaml docstring for the conversion
    from KU-3.17's 100 Pa·s·μm label). We also report the literal
    bead-Stokes value 6 π η_water · bead_radius for diagnostic
    comparison only — it is **not** used in the dynamics.

    Stretch relaxation:        τ_stretch = γ_b ℓ₀ / μ                  [s]
    Bending relaxation:        τ_bend    = γ_b ℓ₀³ / κ                 [s]
    XL relaxation:             τ_xl      = γ_b / k_xl                  [s]
    Tension relaxation:        τ_tension = γ_b R² / (γ_cortex ℓ₀)      [s]
        — derived from f_tension ≈ γ_cortex ℓ₀ / R per bead and the
        local stiffness |∂f/∂R| = γ_cortex ℓ₀ / R², so τ = γ_b/|∂f/∂R|.
    CFL bound:                 dt < α · τ_min, α = cfl_safety_factor.

Free-space convention: a single isolated cell has no periodic box.
We supply ``box_size = 100 · R_cell`` to the ECM mechanics kernels so
the minimum-image arithmetic is a no-op for any realistic cortex bond
or cross-link displacement (≪ box_size / 2).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import yaml


def resolve_cell(cfg: dict[str, Any]) -> dict[str, Any]:
    """Fill the ``cell.derived.*`` block of a Phase 1 Unit 3 config in place."""
    cell = cfg["cell"]
    R = float(cell["R_cell"])
    L_f = float(cell["L_cortex_fiber"])
    n_f = int(cell["n_cortex_fibers"])
    beads = int(cell["beads_per_fiber"])

    cost_cap = int(cell["acceptance"]["n_fibers_max"])
    if n_f > cost_cap and not bool(cell.get("demo_mode", False)):
        raise ValueError(
            f"n_cortex_fibers={n_f} exceeds cost cap {cost_cap}. "
            f"Set cell.demo_mode=true to override."
        )

    circumference = 2.0 * math.pi * R
    coverage_ratio = n_f * L_f / circumference
    rest_length = L_f / (beads - 1)
    backbone_z = (2.0 * (beads - 2) + 2.0) / beads if beads >= 2 else 0.0

    eta = float(cell["water_viscosity"])
    a = float(cell["bead_radius"])
    mu = float(cell["stretching_modulus"])
    kappa = float(cell["bending_modulus"])
    k_xl = float(cell["xl_stiffness"])
    gamma_cortex = float(cell["gamma_cortex"])
    gamma_b = float(cell["gamma_drag_per_bead"])

    gamma_b_water_diag = 6.0 * math.pi * eta * a   # diagnostic only
    tau_xl = gamma_b / k_xl
    tau_stretch = gamma_b * rest_length / mu
    tau_bend = gamma_b * (rest_length ** 3) / kappa
    tau_tension = gamma_b * (R ** 2) / (gamma_cortex * rest_length)
    tau_min = min(tau_xl, tau_stretch, tau_bend, tau_tension)

    derived = cell.setdefault("derived", {})
    derived["cell_circumference"] = circumference
    derived["fiber_coverage_ratio"] = coverage_ratio
    derived["rest_length"] = rest_length
    derived["backbone_z"] = backbone_z
    derived["gamma_b"] = gamma_b
    derived["gamma_b_water_diag"] = gamma_b_water_diag
    derived["tau_xl"] = tau_xl
    derived["tau_stretch"] = tau_stretch
    derived["tau_bend"] = tau_bend
    derived["tau_tension"] = tau_tension
    derived["tau_min"] = tau_min
    derived["box_size"] = 100.0 * R

    dyn = cell.setdefault("dynamics", {})
    alpha = float(dyn.get("cfl_safety_factor", 0.1))
    dt_cfl = alpha * tau_min
    if dyn.get("dt") is None:
        dyn["dt"] = 0.99 * dt_cfl
    derived["dt_cfl"] = dt_cfl

    return cfg


def load_cell_config(path: str | Path) -> dict[str, Any]:
    with open(path, "r") as f:
        cfg = yaml.safe_load(f)
    return resolve_cell(cfg)
