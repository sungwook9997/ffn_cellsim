"""Dynamic osmotic / volume regulation — time-dependent setpoint on the
EXISTING enclosed-volume turgor force (RVD / RVI water-flux kinetics).

ADDITIVE, DEFAULT-OFF module. This module adds **no new force, no new
particles, and no new bonds**. It is a batch :class:`hoomd.custom.Action`
("Updater") that, each batch tick, advances the *setpoint* of the already-built
:class:`ffn_sim.cortex.enclosed_volume.EnclosedVolumePressure` compartment —
i.e. it mutates ``force.p.V0(t)`` and/or ``force.p.turgor_dP0(t)`` in place
along a cited water-flux law. The enclosed-volume force then picks up the new
setpoint on its next ``set_forces`` call. The physics being modelled is
**ion-pump + aquaporin water flux** — regulatory volume decrease/increase
(RVD / RVI): when an osmotic perturbation swells or shrinks the cell, the
plasma membrane's finite water permeability ``Lp`` relaxes the cell volume
back toward an osmotic equilibrium over seconds-to-minutes (Hoffmann 2009).

Why a setpoint Updater and not a lumped force
---------------------------------------------
The HARD full-fidelity rule forbids a lumped proxy force that *replaces*
explicit particles. The mechanical force on every cortex bead remains the
explicit per-bead pressure force of ``EnclosedVolumePressure`` — an explicit
particle/force already in the system. What changes in time is its osmotic
*boundary condition* (the reference volume ``V0`` the cytoplasm bulk law
relaxes toward, and the resting turgor ``Π₀``). Modulating a boundary
condition along a measured transport law (Kedem-Katchalsky / van't Hoff water
flux) is the fine-grained representation of aquaporin/ion-pump regulation at
this ×40 mesoscale — there is no coarser "volume spring" standing in for the
shell. The shell beads still carry the load; we only move the target they
relax to, exactly as a real cell's solute content moves its osmotic setpoint.

Mechanism (explicit transport law)
----------------------------------
The membrane water flux follows the Kedem-Katchalsky / van't Hoff form. The
trans-membrane driving pressure is the difference between the osmotic pressure
the solutes *want* (the regulated target ``ΔP_target``) and the hydrostatic
pressure the cortex currently exerts (the enclosed-volume force's ``ΔP(V)``)::

    Π_osm = R_gas · T · Δc            van't Hoff (Δc = osmolar imbalance)
    dV/dt = -Lp · A_mem · (ΔP_mech(V) − ΔP_target)

with ``Lp`` the membrane hydraulic permeability [m/(s·Pa)] (Olbrich 2000) and
``A_mem`` the membrane area [m²]. Water leaves when the cortex pushes harder
than the osmotic target (cell shrinks toward equilibrium → RVD) and enters
when it pushes less (cell swells → RVI). Discretised over a batch interval
``Δt_batch = batch_steps · dt`` (explicit forward Euler on the slow water
mode)::

    V0(t+Δt) = V0(t) − Lp · A_mem · (ΔP_mech − ΔP_target) · Δt_batch

The relaxation is a first-order exponential toward the target volume. The
restoring stiffness of the WATER mode is the **osmotic modulus** ``Π_osm =
R_gas·T·c_phys`` (for van't Hoff ``Π = R·T·N/V``, ``|V·dΠ/dV| = Π``) — NOT the
cortex mechanical bulk modulus ``K_vol``. The time constant is therefore
``τ_RVD = V0 / (Lp · A_mem · Π_osm)`` (derived in the Sanity Gate), which lands
in the literature seconds-to-minutes RVD/RVI band (Hoffmann 2009) for the
Olbrich ``Lp`` and physiological intracellular osmolarity ``c_phys ≈ 300
mOsm/L`` — a falsifiable consistency check, not a fit.

This module exposes the OSMOTIC TARGET ``ΔP_target`` (or equivalently a target
volume ``V_target``) and the transport coefficients; the Updater integrates
the water mode and writes the new ``V0`` (and optionally a time-varying
``turgor_dP0``) back onto the live ``EnclosedVolumePressure`` force. Setting a
perturbation (e.g. a hyperosmotic shock ``Δc < 0``) and watching ``V0`` relax
is the RVD/RVI experiment.

Sanity Gate
-----------
*Per CLAUDE.md Hard Rule. STATIC + analytical checks in
``ffn_sim/tests/test_osmotic_regulation.py``.*

1. **Dimensional analysis**
   - van't Hoff: ``Π = R_gas·T·Δc`` → [J/(mol·K)]·[K]·[mol/m³] = [J/m³ = Pa].✓
   - Flux: ``Lp`` [m/(s·Pa)]; ``A_mem`` [m²]; ``ΔP`` [Pa] →
     ``dV/dt = Lp·A·ΔP`` = [m/(s·Pa)]·[m²]·[Pa] = [m³/s].✓
   - ``Δt_batch = batch_steps·dt`` [s]; ``dV/dt · Δt`` = [m³].✓
   - ``Π_osm = R_gas·T·c_phys`` = [J/(mol·K)]·[K]·[mol/m³] = [Pa].✓
   - ``τ_RVD = V0/(Lp·A·Π_osm)`` = [m³] / ([m/(s·Pa)]·[m²]·[Pa]) = [s].✓

2. **Boundary cases**
   - ``Lp ≤ 0``, ``A_mem ≤ 0``, ``R_cell ≤ 0``, ``temperature_K ≤ 0``:
     ValueError at resolve.
   - ``enabled: false`` (or missing): resolver returns
     ``ResolvedOsmoticRegulation(enabled=False, …)`` with zeroed transport
     coefficients; the builder returns ``(None, None)`` and attaches nothing
     (OFF-identity — the enclosed-volume setpoint never moves).
   - Already at target (``ΔP_mech == ΔP_target``): ``dV/dt = 0`` — fixed point,
     no setpoint drift.
   - ``Δc`` unset → the osmotic target defaults to the resting turgor already
     resolved on the enclosed-volume force (no spurious perturbation).
   - A floor ``V0_min > 0`` clamps the volume so the shell can never be driven
     to/through zero (avoids the ``R_mean→0`` singular branch in the
     enclosed-volume estimator).

3. **Conservation / no-net-force**
   - This Updater injects NO force and moves NO particle. It only rewrites a
     scalar setpoint on an existing force, so it cannot inject net momentum or
     break Newton's third law. The enclosed-volume force keeps its own
     ``ΣF≈0`` centroid-normal property unchanged (Sanity Gate §3 of
     ``enclosed_volume.py``). STATIC: the Updater's ``act`` leaves
     ``snap.particles.position`` untouched.

4. **Numerical sanity (CFL / stiff note)**
   - The water mode is a SLOW relaxation: ``τ_RVD`` is seconds-to-minutes
     (Hoffmann 2009), i.e. ~10⁸–10¹⁰ × the BAOAB ``dt`` (~10⁻⁸ s). The
     explicit-Euler stability bound on the setpoint mode is
     ``Δt_batch ≤ 2·τ_RVD``; with any sane batch interval this is satisfied
     by ~9 orders of magnitude, so the setpoint integration is
     unconditionally stable here and does NOT tighten the mechanical CFL.
     The Updater asserts ``Δt_batch < τ_RVD`` (well inside the bound) at
     construction and raises otherwise.
   - Forward Euler on V0 is first-order accurate; the slow timescale makes
     the per-batch step error negligible (``Δt_batch/τ_RVD ≪ 1``).

5. **Sign / sense**
   - Cortex pushes harder than osmotic target (``ΔP_mech > ΔP_target``):
     water LEAVES, ``dV/dt < 0`` → ``V0`` decreases → the cell shrinks toward
     equilibrium (RVD after hypotonic swelling). STATIC.
   - Cortex pushes less than the osmotic target (``ΔP_mech < ΔP_target``):
     water ENTERS, ``dV/dt > 0`` → ``V0`` increases → the cell swells toward
     equilibrium (RVI after hypertonic shrinkage). STATIC.
   - ``ΔP_mech == ΔP_target``: ``dV/dt = 0`` (fixed point). STATIC.

6. **Measurement-protocol consistency**
   - The relaxation half-time read from a ``V0(t)`` trajectory must equal
     ``τ_RVD·ln2`` for the chosen ``Lp`` — and ``τ_RVD`` must fall in the
     cited seconds-to-minutes RVD/RVI window. ANALYTICAL test (no long sim):
     one explicit-Euler step on the resolver's transport law reproduces the
     closed-form ``ΔV`` to float precision.

Compartment Performance Contract
--------------------------------
* Particle types added: NONE (modulates an existing force setpoint).
* Particle count at n_fil=1000 and native ~38000 scale: +0 (no particles).
* Bond / angle count: +0 (no bonds, no angles → GAMMA_DENYLIST_PREFIX = '').
* Per-step force: NO. (No force compute of its own; the existing
  ``EnclosedVolumePressure`` remains the only per-step force.)
* Per-batch updater: YES — a ``hoomd.custom.Action`` on a
  ``hoomd.trigger.Periodic(batch_steps)`` trigger. O(1) per tick: it reads
  the force's last-computed ``ΔP``/``V`` diagnostics (already cached by
  ``set_forces``) and writes two scalars. No per-particle loop.
* Uses cpu_local_snapshot: NO (reads the force's cached scalar diagnostics;
  optionally one ``get_snapshot`` only if asked to re-measure V — off by
  default). No GPU↔CPU sync in the hot path.
* Uses cKDTree / broad-phase: NO.
* Hot-path priority: P2 (slow regulatory mode, batched; not in the BAOAB
  inner loop).
* GPU path now: builtin (CPU-side scalar update; nothing to port — it touches
  no per-particle array on the hot path).
* Native ForceCompute candidate: NO (it is not a force; it is a boundary-
  condition scheduler). The force it modulates already has a native/CuPy
  candidate of its own (``enclosed_volume_gpu.py``); this Updater is agnostic
  to which backend computes the pressure.
* Bottleneck risk: negligible — two scalar writes per batch tick.

References
----------
- van't Hoff osmotic law ``Π = R_gas·T·Δc``: van't Hoff 1887; standard
  physical chemistry. ``R_gas = 8.314 462 618 J/(mol·K)`` (CODATA 2018).
- Membrane water permeability ``Lp ~ 1e-12 – 1e-13 m/(s·Pa)``: Olbrich,
  Rawicz, Needham & Evans 2000 (Biophys. J. 79:321, "Water permeability and
  mechanical strength of polyunsaturated lipid bilayers") report bilayer
  ``Pf`` → ``Lp`` of order 1e-12 m/(s·Pa); cells with aquaporins fall in the
  1e-12–1e-13 m/(s·Pa) range (Dvorak; Verkman aquaporin reviews).
- RVD / RVI timescales (seconds-to-minutes): Hoffmann, Lambert & Pedersen
  2009 (Physiol. Rev. 89:193, "Physiology of Cell Volume Regulation in
  Vertebrates") — regulatory volume decrease/increase relax over
  seconds to minutes.
- Intracellular osmolarity ``c_phys ≈ 290–300 mOsm/L`` (= 290–300 mol/m³):
  Lodish et al., Molecular Cell Biology; Alberts et al., Molecular Biology of
  the Cell — standard mammalian intracellular/extracellular osmolarity. This
  sets the osmotic modulus ``Π_osm = R_gas·T·c_phys ≈ 7.7e5 Pa`` (310 K).
- The force being modulated: ``ffn_sim/cortex/enclosed_volume.py``
  (EnclosedVolumePressure; ``V0``, ``turgor_dP0``, ``K_vol`` setpoints).
- KU-3.1 osmotic / turgor anchor: Stewart et al. 2011 (Nature) ΔP ≈ 40 Pa
  interphase → ~400 Pa metaphase (the resting turgor the regulation relaxes
  toward by default).
- Structural analog (batch ``hoomd.custom.Action`` + ``CustomUpdater``):
  ``ffn_sim/cortex/crosslinkers.py`` (XlinkBondUpdater) and
  ``ffn_sim/cell/cell.py`` (SubstrateLigandPin).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import hoomd

from ffn_sim.cortex.enclosed_volume import (
    EnclosedVolumePressure,
    ResolvedEnclosedVolume,
)

# This compartment adds NO bonds (it modulates an existing force setpoint), so
# there is no bond-type prefix to denylist from the cortical-tension estimator.
GAMMA_DENYLIST_PREFIX: str = ""

# CODATA 2018 universal gas constant [J/(mol·K)] — used only by the optional
# van't Hoff target helper; the runtime transport law uses pressures directly.
R_GAS: float = 8.314462618  # J/(mol·K)

# Open PI decisions for this module (empty = none outstanding).
PI_DECISIONS: list[str] = [
    # Lp is given a literature default band [1e-13, 1e-12] m/(s·Pa) (Olbrich "
    # 2000; aquaporin reviews) but is NOT cell-type-specific for MCF7. If a "
    # measured MCF7 plasma-membrane Lp is required for production, surface to "
    # PI rather than picking a value inside the band.
    "Lp default = 1e-13 m/(s·Pa) (Olbrich 2000 band low end, aquaporin-poor "
    "membrane); no MCF7-specific Lp datum found — PI to confirm or supply a "
    "measured value before production RVD/RVI runs.",
]


# ---------------------------------------------------------------------------
# Resolved parameters
# ---------------------------------------------------------------------------
@dataclass(slots=True)
class ResolvedOsmoticRegulation:
    """Dynamic osmotic / volume-regulation parameters (all SI).

    Off-identity: ``enabled=False`` => every transport coefficient is zero,
    the builder attaches nothing, and the enclosed-volume setpoint never moves.
    """

    enabled: bool                       # —    master switch (default-OFF)

    # Transport law (Kedem-Katchalsky / van't Hoff water flux).
    Lp: float                           # m/(s·Pa)  membrane hydraulic permeability
    A_mem: float                        # m²        membrane area (water-exchange area)
    batch_steps: int                    # —    integration steps between setpoint updates
    dt: float                           # s    host BAOAB timestep
    batch_dt: float                     # s    Δt_batch = batch_steps · dt (derived)

    # Regulated osmotic target. The setpoint relaxes the enclosed-volume force
    # toward the pressure ``dP_target`` (Pa). Defaults to the resting turgor of
    # the enclosed-volume force (no perturbation) when no Δc is supplied.
    dP_target: float                    # Pa   osmotic target pressure (van't Hoff or turgor)
    temperature_K: float                # K    absolute temperature (van't Hoff)

    # Safety floor on the regulated reference volume (Pa-side singular guard).
    V0_min: float                       # m³   minimum allowed V0 (R_mean→0 guard)

    # Derived diagnostics.
    K_vol: float = 0.0                  # Pa   cortex bulk modulus copied from the EV force
    Pi_osm: float = 0.0                 # Pa   osmotic modulus = R_gas·T·c_phys (restoring stiffness)
    V0_ref: float = 0.0                 # m³   initial reference volume (relax start)
    tau_RVD: float = 0.0                # s    relaxation time V0/(Lp·A·Pi_osm)

    # Physiological intracellular osmolarity (van't Hoff modulus anchor).
    c_phys: float = 0.0                 # mol/m³  total intracellular osmolyte conc.

    # Optional explicit osmolar imbalance (van't Hoff source of dP_target).
    delta_c: float | None = None        # mol/m³  Δ solute concentration (None=unset)


def _require_finite_positive(name: str, x: float) -> None:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} must be finite and > 0; got {x!r}")


def vant_hoff_pressure(delta_c: float, temperature_K: float) -> float:
    """van't Hoff osmotic pressure ``Π = R_gas·T·Δc`` [Pa].

    Acceptance helper. ``delta_c`` is the osmolar imbalance [mol/m³] and
    ``temperature_K`` the absolute temperature [K].

    Args:
        delta_c: Solute concentration imbalance across the membrane [mol/m³].
            Positive => more solute inside => water enters => swelling driver.
        temperature_K: Absolute temperature [K].

    Returns:
        Osmotic pressure [Pa].
    """
    return R_GAS * temperature_K * delta_c


def resolve_osmotic_regulation(
    cfg: dict,
    *,
    R_cell: float,
    dt: float,
    p_enclosed_volume: ResolvedEnclosedVolume | None = None,
    temperature_K: float = 310.15,
) -> ResolvedOsmoticRegulation:
    """Resolve the ``osmotic_regulation`` config block.

    Default-OFF: a missing block, or ``enabled: false``, returns a disabled
    resolved object with zeroed transport coefficients (the builder then
    attaches nothing and the enclosed-volume setpoint never moves).

    The transport law and timescale ``τ_RVD = V0/(Lp·A_mem·Π_osm)`` use the
    OSMOTIC modulus ``Π_osm = R_gas·T·c_phys`` (the restoring stiffness of the
    WATER mode) — NOT the cortex mechanical bulk modulus ``K_vol`` (the latter
    is copied from the host force only as a diagnostic). The host
    enclosed-volume force's ``V0`` (passed via ``p_enclosed_volume``) sets the
    relaxation start so the regulation is consistent with the compartment it
    modulates. With no ``p_enclosed_volume`` the timescale is left at 0
    (diagnostic-only) and the membrane area defaults to the cell's sphere area
    ``4π R_cell²``.

    Args:
        cfg: YAML root, the ``cortex:`` sub-dict, or the
            ``cortex.osmotic_regulation`` sub-dict.
        R_cell: Nominal cell radius [m] (sets the default membrane area and the
            default reference volume).
        dt: Host BAOAB timestep [s] (for ``batch_dt`` and the slow-mode CFL).
        p_enclosed_volume: The resolved enclosed-volume force parameters whose
            setpoint this module modulates. Required to compute ``τ_RVD`` and
            to default the osmotic target to the resting turgor.
        temperature_K: Absolute temperature [K] (default 310.15 K = 37 °C,
            physiological — van't Hoff factor). Overridable via cfg.

    Returns:
        ResolvedOsmoticRegulation.

    Raises:
        ValueError: on negative/zero transport coefficients, bad batch_steps,
            or a slow-mode CFL violation.
    """
    if "cortex" in cfg:
        cfg = cfg["cortex"]
    if "osmotic_regulation" in cfg:
        cfg = cfg["osmotic_regulation"]

    enabled = bool(cfg.get("enabled", False))

    # Defaults sourced from the host enclosed-volume force when present.
    K_vol = float(p_enclosed_volume.K_vol) if p_enclosed_volume is not None else 0.0
    V0_ref = (
        float(p_enclosed_volume.V0)
        if p_enclosed_volume is not None
        else (4.0 / 3.0) * math.pi * R_cell**3
    )
    turgor0 = (
        float(p_enclosed_volume.turgor_dP0) if p_enclosed_volume is not None else 0.0
    )

    if not enabled:
        # OFF-identity: zeroed transport coefficients, nothing to attach.
        return ResolvedOsmoticRegulation(
            enabled=False,
            Lp=0.0,
            A_mem=0.0,
            batch_steps=0,
            dt=float(dt),
            batch_dt=0.0,
            dP_target=0.0,
            temperature_K=float(cfg.get("temperature_K", temperature_K)),
            V0_min=0.0,
            K_vol=K_vol,
            Pi_osm=0.0,
            V0_ref=V0_ref,
            tau_RVD=0.0,
            c_phys=0.0,
            delta_c=None,
        )

    # --- enabled: validate the physical transport parameters ---------------
    _require_finite_positive("R_cell", R_cell)
    _require_finite_positive("dt", float(dt))

    temperature_K = float(cfg.get("temperature_K", temperature_K))
    _require_finite_positive("temperature_K", temperature_K)

    # Membrane hydraulic permeability. Literature band [1e-13, 1e-12]
    # m/(s·Pa) (Olbrich 2000; aquaporin reviews). Default = low end (PI_DECISIONS).
    Lp = float(cfg.get("Lp", 1.0e-13))  # m/(s·Pa)  Olbrich 2000 band low end
    _require_finite_positive("Lp", Lp)

    # Membrane / water-exchange area; defaults to the cell sphere area.
    A_mem = float(cfg.get("A_mem", 4.0 * math.pi * R_cell**2))  # m²
    _require_finite_positive("A_mem", A_mem)

    batch_steps = int(cfg.get("batch_steps", 1000))
    if batch_steps <= 0:
        raise ValueError(f"batch_steps must be a positive int; got {batch_steps}")
    batch_dt = batch_steps * float(dt)

    # Osmotic target. Either an explicit Δc (van't Hoff) → ΔP_target relative to
    # the resting turgor, or an explicit dP_target, or default = resting turgor
    # (no perturbation; fixed point at construction volume).
    delta_c = cfg.get("delta_c", None)
    if delta_c is not None:
        delta_c = float(delta_c)
        if not math.isfinite(delta_c):
            raise ValueError(f"delta_c must be finite; got {delta_c!r}")
        # van't Hoff perturbation pressure added to the resting turgor target.
        dP_target = turgor0 + vant_hoff_pressure(delta_c, temperature_K)
    else:
        dP_target = float(cfg.get("dP_target", turgor0))
    if not math.isfinite(dP_target):
        raise ValueError(f"dP_target must be finite; got {dP_target!r}")

    # Volume floor (singular-guard for the R_mean→0 estimator branch). Default
    # to a small fraction of the reference volume (grid-invariant; not tuned).
    V0_min = float(cfg.get("V0_min", 0.01 * V0_ref))
    _require_finite_positive("V0_min", V0_min)
    if V0_min >= V0_ref:
        raise ValueError(
            f"V0_min ({V0_min:.3e} m³) must be < V0_ref ({V0_ref:.3e} m³)."
        )

    # Physiological intracellular osmolarity (van't Hoff modulus anchor). The
    # OSMOTIC restoring stiffness of the water mode is the osmotic modulus
    # Π_osm = R_gas·T·c_phys (for van't Hoff Π = R·T·N/V, |V·dΠ/dV| = Π) — NOT
    # the cortex mechanical bulk modulus K_vol. Default c_phys = 300 mol/m³
    # (= 300 mOsm/L), the standard mammalian intracellular osmolarity
    # (Lodish, Molecular Cell Biology; Alberts MBoC: ~290-300 mOsm). Grid-
    # invariant intensive concentration, not a fit.
    c_phys = float(cfg.get("c_phys", 300.0))  # mol/m³  intracellular osmolarity
    _require_finite_positive("c_phys", c_phys)
    Pi_osm = vant_hoff_pressure(c_phys, temperature_K)  # Pa  osmotic modulus

    # Relaxation timescale τ_RVD = V0 / (Lp·A_mem·Π_osm). The osmotic modulus
    # Π_osm (~7.7e5 Pa at 300 mOsm, 310 K) is the restoring stiffness of the
    # water mode, giving the Hoffmann 2009 seconds-to-minutes RVD/RVI window.
    tau_RVD = V0_ref / (Lp * A_mem * Pi_osm)
    # Slow-mode explicit-Euler stability: Δt_batch ≤ 2·τ_RVD (we demand the
    # stronger Δt_batch < τ_RVD so the regulation is well-resolved).
    if batch_dt >= tau_RVD:
        raise ValueError(
            f"osmotic_regulation slow-mode CFL violated: batch_dt = "
            f"{batch_dt:.3e} s ≥ τ_RVD = {tau_RVD:.3e} s "
            f"(τ_RVD = V0/(Lp·A·Π_osm), Lp={Lp:.2e} m/(s·Pa), "
            f"A_mem={A_mem:.3e} m², Π_osm={Pi_osm:.3e} Pa). "
            "Reduce batch_steps (resolve the slow water mode more finely) "
            "or revisit Lp/A_mem against the literature band."
        )

    return ResolvedOsmoticRegulation(
        enabled=True,
        Lp=Lp,
        A_mem=A_mem,
        batch_steps=batch_steps,
        dt=float(dt),
        batch_dt=batch_dt,
        dP_target=dP_target,
        temperature_K=temperature_K,
        V0_min=V0_min,
        K_vol=K_vol,
        Pi_osm=Pi_osm,
        V0_ref=V0_ref,
        tau_RVD=tau_RVD,
        c_phys=c_phys,
        delta_c=delta_c,
    )


# ---------------------------------------------------------------------------
# Closed-form transport oracle (used by the Updater AND by tests)
# ---------------------------------------------------------------------------
def water_flux_volume_step(
    p: ResolvedOsmoticRegulation,
    *,
    dP_mech: float,
    V0: float,
) -> float:
    """One explicit-Euler step of the membrane water-flux law on ``V0``.

    Implements ``V0(t+Δt) = V0 − Lp·A_mem·(ΔP_mech − ΔP_target)·Δt_batch`` and
    clamps to the volume floor ``V0_min``.

    Sign sense: when the cortex pushes harder than the osmotic target
    (``dP_mech > dP_target``) water LEAVES and ``V0`` decreases (RVD); when it
    pushes less, water ENTERS and ``V0`` increases (RVI).

    Args:
        p: Resolved osmotic-regulation parameters.
        dP_mech: Current mechanical (enclosed-volume) pressure ΔP(V) [Pa],
            i.e. the hydrostatic pressure the cortex currently exerts.
        V0: Current reference volume [m³].

    Returns:
        Updated reference volume [m³] (clamped at ``V0_min``).
    """
    dVdt = -p.Lp * p.A_mem * (dP_mech - p.dP_target)  # m³/s
    V0_new = V0 + dVdt * p.batch_dt
    return max(V0_new, p.V0_min)


# ---------------------------------------------------------------------------
# Batch Updater Action: advance the enclosed-volume setpoint over time
# ---------------------------------------------------------------------------
class OsmoticRegulationUpdater(hoomd.custom.Action):
    """Batch Action that steps the enclosed-volume setpoint along the
    water-flux law (RVD / RVI).

    Each tick (``hoomd.trigger.Periodic(batch_steps)``):

    1. Read the live mechanical pressure ``ΔP_mech`` the enclosed-volume force
       last computed (its cached ``last_pressure`` diagnostic — already in
       float64, no extra snapshot/sync). On the first tick before any
       ``set_forces`` has run, fall back to the resting turgor.
    2. Take one explicit-Euler water-flux step on ``V0``
       (:func:`water_flux_volume_step`).
    3. Write the new ``V0`` back onto the live force's ``.p`` dataclass
       (mutable, slots) so the next ``set_forces`` relaxes the shell toward the
       updated reference volume. NO particle is moved, NO force is added.

    The Action mutates only two scalars on the force; it touches no
    per-particle array, so it injects no net momentum and never breaks
    Newton's third law (Sanity Gate §3).

    Parameters
    ----------
    p : ResolvedOsmoticRegulation
        Resolved regulation parameters (must be ``enabled``).
    ev_force : EnclosedVolumePressure
        The live enclosed-volume force whose ``V0`` setpoint is modulated.
    """

    def __init__(
        self,
        p: ResolvedOsmoticRegulation,
        ev_force: EnclosedVolumePressure,
    ) -> None:
        super().__init__()
        if not p.enabled:
            # Defensive: a disabled config should never reach the Updater (the
            # builder early-returns), but make the no-op explicit.
            raise ValueError(
                "OsmoticRegulationUpdater constructed with a disabled config; "
                "the builder must early-return for enabled=False."
            )
        self.p = p
        self.ev_force = ev_force
        self._sim_ref: hoomd.Simulation | None = None
        # Diagnostics for tests / monitors.
        self.last_V0: float = float(ev_force.p.V0)
        self.last_dP_mech: float = float("nan")
        self.n_ticks: int = 0

    def attach(self, simulation: hoomd.Simulation) -> None:  # noqa: D401
        super().attach(simulation)
        self._sim_ref = simulation

    def act(self, timestep: int) -> None:  # noqa: D401
        # Current mechanical pressure from the force's cached diagnostic. If
        # the force has not yet run a set_forces (last_pressure is NaN), use the
        # resting turgor as the starting ΔP so the first step is well-defined.
        dP_mech = float(self.ev_force.last_pressure)
        if not math.isfinite(dP_mech):
            dP_mech = float(self.ev_force.p.turgor_dP0)

        V0_old = float(self.ev_force.p.V0)
        V0_new = water_flux_volume_step(self.p, dP_mech=dP_mech, V0=V0_old)

        # Write the new setpoint back onto the live force (mutable dataclass).
        self.ev_force.p.V0 = V0_new

        self.last_dP_mech = dP_mech
        self.last_V0 = V0_new
        self.n_ticks += 1


# ---------------------------------------------------------------------------
# Public builder: attach the regulation Updater to a built simulation
# ---------------------------------------------------------------------------
def attach_osmotic_regulation_to_simulation(
    sim: hoomd.Simulation,
    p: ResolvedOsmoticRegulation,
    ev_force: EnclosedVolumePressure | None,
) -> tuple[OsmoticRegulationUpdater | None, hoomd.update.CustomUpdater | None]:
    """Attach the regulation Updater to an already-built simulation.

    OFF-identity: when ``p.enabled`` is False (or ``ev_force`` is None) this
    returns ``(None, None)`` and attaches NOTHING — the simulation/snapshot is
    unchanged and the enclosed-volume setpoint never moves.

    Args:
        sim: Built simulation that already carries the enclosed-volume force.
        p: Resolved osmotic-regulation parameters.
        ev_force: The live :class:`EnclosedVolumePressure` force whose ``V0``
            setpoint will be modulated. If None, nothing is attached.

    Returns:
        ``(action, custom_updater)`` when enabled, else ``(None, None)``.
    """
    if not p.enabled or ev_force is None:
        return None, None

    action = OsmoticRegulationUpdater(p, ev_force)
    updater = hoomd.update.CustomUpdater(
        action=action, trigger=hoomd.trigger.Periodic(p.batch_steps)
    )
    sim.operations.updaters.append(updater)
    return action, updater
