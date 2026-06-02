"""Resolved parameter set for the Layer-2 CBM spheroid (HOOMD-free; testable standalone).

``resolve_layer2(cfg)`` reads the primary, literature-anchored scales from a config dict
(the ``spheroid:`` block of ``configs/layer2_cbm.yaml``) and computes every derived
runtime quantity (cell radius, Morse parameters, per-cell Stokes drag, CFL timestep) so
the HOOMD builder in ``cbm.py`` consumes only resolved SI numbers. Pure python — no HOOMD,
no I/O — so the derivations are unit-testable without a simulation.

Magic-Number Block (CLAUDE.md hard rule)
----------------------------------------
Every primary is literature-anchored or a flagged modeling choice; every derived value is
COMPUTED from primaries (never hand-set):

- ``diameter`` 15 µm  -- MCF7 suspended single-cell, Wagner 2011 (Coulter, MEASURED,
  PMC3147247: 14.8 µm); 15 µm = low end of the measured 15-20 µm band (band-position is a
  documented modeling choice). See ``docs/LAYER2_ANCHORS_2026-06-02.md``.
- ``D_e = 2 * F_detach * contact_zone_width`` -- cell-cell adhesion well depth, set so the
  Morse MAX attractive force (D_e*alpha/2) equals the MEASURED MCF7-MCF7 de-adhesion force.
  Source: Iturri 2020 *Cells* PMC7227807 (MCF7-MCF7 SCFS, MEASURED, open-access): detachment
  force ~6-7 nN at 120 s mature contact (~1-2 nN at 5 s nascent); corroborated by Omidvar
  2014 *J Biomech* (ranking MCF7>T47D>MDA, nN-scale) + Omidvar 2016 *Mol Cell Biochem*
  (1.42-2.85 nN). This is the nN-scale quantity that balances nN single-cell traction and
  sets when cells separate during spreading. NOTE: this RETIRES the earlier cadherin
  ``N_cad*<F>*x_beta`` construction (=1.2e-17 J, ~4 orders of magnitude too small — a pN
  single-bond product) and the misattributed Buckley-2014 Δx*=4nm (that is the αE-catenin–
  actin bond, NOT E-cadherin homophilic; the real E-cadherin catch-bond is Rakshit 2012
  PNAS, deferred to the L2.5 KU-4.2 upgrade). The measured work-of-de-adhesion (~200 fJ)
  is a scale cross-check; it exceeds 2*F_detach*contact_zone because SCFS work integrates
  membrane-tether pulling beyond the CBM contact range (~µm).
- ``gamma_cell = 6*pi*eta*R_cell`` -- per-cell Stokes drag (KU-1.26), DERIVED; matches the
  single-cell BAOAB drag form (``cell/dt_reconcile.py``).
- ``morse_r0 = diameter`` -- two R_cell spheres touch (surfaces) at center-to-center 2R.
- ``morse_alpha = 1/contact_zone_width`` -- inverse adhesive range; contact_zone_width is a
  geometric modeling choice (~10% of diameter), not a literature constant.
- ``dt_cfl = safety * gamma/k_spring`` with ``k_spring = 2*D_e*alpha^2`` (Morse curvature
  at the minimum) -- overdamped relaxation-time CFL.

Sanity Gate
-----------
- Dimensional: diameter/widths [m]; force [N]; x_beta [m] -> D_e [J]; eta [Pa·s] ->
  gamma [N·s/m]; dt [s]. Verified by construction + tests.
- All primaries and derived values are finite and strictly positive (guarded).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = ["ResolvedL2", "resolve_layer2"]


@dataclass(frozen=True)
class ResolvedL2:
    """Fully-resolved Layer-2 CBM parameters (all SI)."""

    # ---- primaries (literature-anchored / flagged modeling choices) ----
    temperature: float          # K
    kT: float                   # J
    water_viscosity: float      # Pa·s
    seed: int
    diameter: float                 # m  (MCF7, Wagner 2011)
    deadhesion_force_mature: float  # N  (measured MCF7-MCF7 120s, Iturri 2020)
    deadhesion_force_nascent: float # N  (measured MCF7-MCF7 5s; maturation range)
    contact_zone_width: float       # m  (geometric modeling choice; Morse range)
    cfl_safety_factor: float        # –

    # ---- derived (computed from primaries) ----
    R_cell: float               # m  = diameter/2
    D_e: float                  # J  = 2·F_detach·contact_zone (Morse max force = measured F_detach)
    morse_r0: float             # m  = diameter
    morse_alpha: float          # 1/m = 1/contact_zone_width
    morse_k_spring: float       # N/m = 2·D_e·alpha²  (curvature at minimum)
    gamma_cell: float           # N·s/m = 6π·η·R_cell
    dt_cfl: float               # s  = safety·gamma/k_spring


def _require_positive(name: str, value: float) -> float:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and strictly positive; got {value!r}.")
    return value


def resolve_layer2(cfg: dict) -> ResolvedL2:
    """Resolve the ``spheroid:`` config block into a fully-derived ``ResolvedL2``.

    Args:
        cfg: parsed YAML dict containing a top-level ``spheroid`` mapping.

    Returns:
        ``ResolvedL2`` with all derived SI quantities computed.

    Raises:
        KeyError: if a required key is missing.
        ValueError: if any primary or derived value is non-finite or non-positive.
    """
    s = cfg["spheroid"]
    cell = s["cell"]
    adh = s["adhesion"]
    dyn = s["dynamics"]

    temperature = _require_positive("temperature", float(s["temperature"]))
    kT = _require_positive("kT", float(s["kT"]))
    eta = _require_positive("water_viscosity", float(s["water_viscosity"]))
    seed = int(s["seed"])

    diameter = _require_positive("cell.diameter", float(cell["diameter"]))
    f_detach = _require_positive(
        "adhesion.deadhesion_force_mature", float(adh["deadhesion_force_mature"])
    )
    f_detach_nascent = _require_positive(
        "adhesion.deadhesion_force_nascent", float(adh["deadhesion_force_nascent"])
    )
    contact_zone = _require_positive(
        "adhesion.contact_zone_width", float(adh["contact_zone_width"])
    )
    safety = _require_positive("dynamics.cfl_safety_factor", float(dyn["cfl_safety_factor"]))

    # ---- derivations ----
    R_cell = diameter / 2.0
    morse_r0 = diameter                              # surfaces touch at center-sep 2R
    morse_alpha = 1.0 / contact_zone                 # inverse adhesive range
    # D_e set so Morse max attractive force (D_e·alpha/2) == measured MCF7-MCF7 de-adhesion
    # force (Iturri 2020) — the nN quantity that balances nN traction. RETIRES the pN-scale
    # N_cad·<F>·x_beta seed (was ~4 OOM too small).
    D_e = 2.0 * f_detach * contact_zone
    morse_k_spring = 2.0 * D_e * morse_alpha * morse_alpha   # Morse curvature at min
    gamma_cell = 6.0 * math.pi * eta * R_cell        # Stokes drag, KU-1.26
    dt_cfl = safety * gamma_cell / morse_k_spring    # overdamped relaxation CFL

    for name, val in [
        ("R_cell", R_cell), ("D_e", D_e), ("morse_r0", morse_r0),
        ("morse_alpha", morse_alpha), ("morse_k_spring", morse_k_spring),
        ("gamma_cell", gamma_cell), ("dt_cfl", dt_cfl),
    ]:
        _require_positive(name, val)

    return ResolvedL2(
        temperature=temperature, kT=kT, water_viscosity=eta, seed=seed,
        diameter=diameter, deadhesion_force_mature=f_detach,
        deadhesion_force_nascent=f_detach_nascent,
        contact_zone_width=contact_zone, cfl_safety_factor=safety,
        R_cell=R_cell, D_e=D_e, morse_r0=morse_r0, morse_alpha=morse_alpha,
        morse_k_spring=morse_k_spring, gamma_cell=gamma_cell, dt_cfl=dt_cfl,
    )
