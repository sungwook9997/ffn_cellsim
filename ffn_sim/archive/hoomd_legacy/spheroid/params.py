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
- ``cortical_tension`` 0.57 mN/m -- physiological single-cell cortical tension input for
  Layer-2 (KU-3.5 g_rigid native, inside the 0.35-0.65 mN/m band). Track A consumes this
  from config only; Track B owns the fine-grained myosin->actin producer.
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
- ``gamma_cell = n_eng * kappa_clutch / k_off`` (n_eng = N_clutch * k_on/(k_on+k_off)) --
  per-cell MIGRATION drag = adhesion-limited motor-clutch friction, DERIVED from the
  single-cell clutch ensemble (scale-bridge from Layer-1; KU-2.18 Bangasser 2013). This
  RETIRES the water-Stokes 6*pi*eta*R (~9.8e-8 N·s/m), which is ~6.5 orders too small for
  crawling (a nN traction would give an unphysical mm/s). Cross-validated against measured
  MCF7 motility (0.2-0.5 um/min under nN traction => ~0.3 N·s/m; Maiuri PMC4245034,
  SICM-TFM PMC8697701). At gamma=0.3, a 1 nN net traction gives v=0.2 um/min (slow epithelial).
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

__all__ = ["ResolvedL2", "resolve_layer2", "ResolvedProliferation", "resolve_proliferation"]


@dataclass(frozen=True)
class ResolvedL2:
    """Fully-resolved Layer-2 CBM parameters (all SI)."""

    # ---- primaries (literature-anchored / flagged modeling choices) ----
    temperature: float          # K
    kT: float                   # J
    water_viscosity: float      # Pa·s
    seed: int
    diameter: float                 # m  (MCF7, Wagner 2011)
    cortical_tension: float         # N/m  Layer-2 input γ (KU-3.5 g_rigid; config, not cortex import)
    cortical_tension_band: tuple[float, float] # N/m  accepted physiological band for γ
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
    gamma_cell: float           # N·s/m  per-cell MIGRATION drag = clutch-ensemble friction
                                #         n_eng·κ_clutch/k_off (KU-2.18 Bangasser 2013); the
                                #         water-Stokes 6π·η·R is RETIRED (see Magic-Number Block)
    dt_cfl: float               # s  = safety·gamma/k_spring


@dataclass(frozen=True)
class ResolvedProliferation:
    """Resolved Layer-2 contact-inhibited proliferation parameters (L2.4).

    Kept SEPARATE from ``ResolvedL2`` so the G1/G3 parameter contract (and its tests) is
    unchanged: proliferation is an additive mechanism layer, resolved on demand.
    """

    cycle_time_mean: float      # s   uncrowded MCF7 doubling (BNID 100685)
    cycle_time_cv: float        # –   per-cell cycle-time CV (modeling choice)
    kissing_number: int         # –   3D sphere kissing number Z=12 (geometric)
    shell_factor: float         # –   first-shell cutoff in units of r0 (geometric)
    split_factor: float         # –   daughter separation in units of r0
    # derived (need r0 from ResolvedL2 → resolved against it)
    shell_cutoff: float         # m   = shell_factor · morse_r0
    split_distance: float       # m   = split_factor · morse_r0
    min_gap: float              # m   free-space threshold = morse_r0 − contact_zone_width
                                #     (Morse repulsive-core onset; the Drasdo-Höhme room test)


def resolve_proliferation(cfg: dict, resolved: "ResolvedL2") -> ResolvedProliferation:
    """Resolve the ``spheroid.proliferation`` block against an already-resolved ``ResolvedL2``.

    Args:
        cfg: parsed YAML dict containing ``spheroid.proliferation``.
        resolved: the resolved CBM parameters (supplies ``morse_r0`` for the length derivations).

    Returns:
        ``ResolvedProliferation`` with the two length scales (shell cutoff, daughter split)
        derived from ``morse_r0``.

    Raises:
        KeyError: if a required key is missing.
        ValueError: if any value is out of its valid range.
    """
    p = cfg["spheroid"]["proliferation"]
    cycle_mean = _require_positive("proliferation.cycle_time_mean", float(p["cycle_time_mean"]))
    cycle_cv = float(p["cycle_time_cv"])
    if not (0.0 <= cycle_cv < 1.0):
        raise ValueError(f"proliferation.cycle_time_cv must be in [0, 1); got {cycle_cv!r}.")
    z = int(p["kissing_number"])
    if z < 1:
        raise ValueError(f"proliferation.kissing_number must be >= 1; got {z!r}.")
    shell_factor = _require_positive("proliferation.shell_factor", float(p["shell_factor"]))
    if shell_factor < 1.0:
        raise ValueError(
            f"proliferation.shell_factor must be >= 1 (>= rest separation); got {shell_factor!r}."
        )
    split_factor = _require_positive("proliferation.split_factor", float(p["split_factor"]))

    # Free-space (room) threshold = Morse repulsive-core onset, one contact-zone inside the
    # rest separation. DERIVED from the pair potential, not tuned. Must stay positive.
    min_gap = _require_positive(
        "proliferation.min_gap (= morse_r0 − contact_zone_width)",
        resolved.morse_r0 - resolved.contact_zone_width,
    )

    return ResolvedProliferation(
        cycle_time_mean=cycle_mean,
        cycle_time_cv=cycle_cv,
        kissing_number=z,
        shell_factor=shell_factor,
        split_factor=split_factor,
        shell_cutoff=shell_factor * resolved.morse_r0,
        split_distance=split_factor * resolved.morse_r0,
        min_gap=min_gap,
    )


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
    mech = s["mechanics"]
    adh = s["adhesion"]
    dyn = s["dynamics"]

    temperature = _require_positive("temperature", float(s["temperature"]))
    kT = _require_positive("kT", float(s["kT"]))
    eta = _require_positive("water_viscosity", float(s["water_viscosity"]))
    seed = int(s["seed"])

    diameter = _require_positive("cell.diameter", float(cell["diameter"]))
    cortical_tension = _require_positive(
        "mechanics.cortical_tension", float(mech["cortical_tension"])
    )
    gamma_band_raw = mech.get("cortical_tension_band", (0.0, float("inf")))
    if len(gamma_band_raw) != 2:
        raise ValueError("mechanics.cortical_tension_band must contain [lo, hi].")
    gamma_band = (
        _require_positive("mechanics.cortical_tension_band[0]", float(gamma_band_raw[0])),
        _require_positive("mechanics.cortical_tension_band[1]", float(gamma_band_raw[1])),
    )
    if gamma_band[0] > gamma_band[1]:
        raise ValueError("mechanics.cortical_tension_band must satisfy lo <= hi.")
    if not (gamma_band[0] <= cortical_tension <= gamma_band[1]):
        raise ValueError(
            "mechanics.cortical_tension must be inside cortical_tension_band; "
            f"got {cortical_tension!r} not in {gamma_band!r}."
        )
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
    clutch_n = _require_positive("dynamics.clutch_n", float(dyn["clutch_n"]))
    clutch_k_on = _require_positive("dynamics.clutch_k_on", float(dyn["clutch_k_on"]))
    clutch_k_off = _require_positive("dynamics.clutch_k_off", float(dyn["clutch_k_off"]))
    clutch_kappa = _require_positive("dynamics.clutch_stiffness", float(dyn["clutch_stiffness"]))

    # ---- derivations ----
    R_cell = diameter / 2.0
    morse_r0 = diameter                              # surfaces touch at center-sep 2R
    morse_alpha = 1.0 / contact_zone                 # inverse adhesive range
    # D_e set so Morse max attractive force (D_e·alpha/2) == measured MCF7-MCF7 de-adhesion
    # force (Iturri 2020) — the nN quantity that balances nN traction. RETIRES the pN-scale
    # N_cad·<F>·x_beta seed (was ~4 OOM too small).
    D_e = 2.0 * f_detach * contact_zone
    morse_k_spring = 2.0 * D_e * morse_alpha * morse_alpha   # Morse curvature at min
    # Per-cell migration drag = motor-clutch ensemble friction (KU-2.18 scale-bridge),
    # NOT water-Stokes (~6.5 OOM too small for crawling). gamma = n_eng * kappa / k_off.
    n_eng = clutch_n * clutch_k_on / (clutch_k_on + clutch_k_off)
    gamma_cell = n_eng * clutch_kappa / clutch_k_off
    dt_cfl = safety * gamma_cell / morse_k_spring    # overdamped relaxation CFL

    for name, val in [
        ("R_cell", R_cell), ("D_e", D_e), ("morse_r0", morse_r0),
        ("morse_alpha", morse_alpha), ("morse_k_spring", morse_k_spring),
        ("gamma_cell", gamma_cell), ("dt_cfl", dt_cfl),
    ]:
        _require_positive(name, val)

    return ResolvedL2(
        temperature=temperature, kT=kT, water_viscosity=eta, seed=seed,
        diameter=diameter, cortical_tension=cortical_tension,
        cortical_tension_band=gamma_band, deadhesion_force_mature=f_detach,
        deadhesion_force_nascent=f_detach_nascent,
        contact_zone_width=contact_zone, cfl_safety_factor=safety,
        R_cell=R_cell, D_e=D_e, morse_r0=morse_r0, morse_alpha=morse_alpha,
        morse_k_spring=morse_k_spring, gamma_cell=gamma_cell, dt_cfl=dt_cfl,
    )
