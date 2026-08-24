"""3-model ERM membrane--cortex connector spectrum for the Active Cell runtime (Warp CUDA source).

This EXTENDS the single-tether primitives in :mod:`aleph.components.incumbent.erm_tether` (``erm_bell_force_kernel`` /
``erm_bell_kmc_kernel`` / ``erm_tether_rupture_kernel``) from a discrete per-tether on/off bond into a
mean-field connector spectrum that brackets the ERM biology.  It keeps the SAME force/commit split those
kernels established: force kernels are pure (an inner quasi-static iteration may discover a loaded connector
but may not advance biology), and the bond-density integrator commits ONLY once at an accepted outer physical
transaction under a device ``accepted`` latch — so a rejected mechanical candidate is a bit-exact no-op.

The spectrum (host algebra + provenance in :mod:`erm_spectrum_analytic`; motivation from Fehon, McClatchey &
Bretscher 2010, *Nat Rev Mol Cell Biol* 11:276 — ERM is a regulated dynamic clutch, not a static weld):

* **model 1 -- RIGID_WELD (완전결합):** ``erm_weld_force_kernel`` — a *bilateral* stiff penalty enforcing the
  Δu = 0 co-motion limit (``k_weld -> inf`` pins the membrane node to its cortex node).
* **model 2 -- ELASTIC_FRICTION (탄성·마찰):** ``erm_clutch_force_kernel`` with a *fixed* bond density —
  a unilateral tension spring ``rho_b·k_b·max(L-rest,0)`` (ezrin bears tension / slides, never pushes) plus a
  viscous along-bond friction ``xi·du_dot`` (Kelvin--Voigt clutch).
* **model 3 -- DYNAMIC_CLUTCH (동적결합):** the SAME ``erm_clutch_force_kernel`` reading an *evolving* bond
  density that ``erm_bond_density_commit_kernel`` advances by the mean-field balance
  ``d(rho_b)/dt = k_on·(rho_max - rho_b) - k_off(f)·rho_b`` on the per-bond load ``f = k_b·max(L-rest,0)``,
  with ``k_off`` either Bell slip (reusing :func:`aleph.components.motor.hand.bell_off_rate`) or Pereverzev
  catch/slip.  This is the continuum generalization of ``erm_bell_kmc_kernel``: the discrete bond is the
  ``rho_max = 1`` single-linker limit.

The kinetic parameters intentionally have NO biological defaults; :class:`ERMSpectrumParams` validates a
complete, source-labelled, per-model contract, mirroring :class:`aleph.components.incumbent.erm_tether.ERMBellKinetics`.

Sanity Gate:
    * Force: model-2/3 springs are unilateral (zero in compression/at rest, tensile Newton-3rd when extended);
      model-1 is bilateral (restores compression too).  A membrane--actin linker is a tether, not a strut,
      except in the explicit rigid-weld limit which is a constraint penalty.
    * Time: force kernels never write ``bound`` or ``rho_b``.  Only the accepted-gated
      ``erm_bond_density_commit_kernel`` advances ``rho_b``; solver iterations are not biological time.
    * Kinetics: catch/slip off-rate falls then rises with per-bond load; ``rho_b`` relaxes to
      ``k_on·rho_max/(k_on+k_off)`` via an exact constant-coefficient tick (unconditionally stable, no
      overshoot); a rejected transaction is an exact no-op.
    * Residency: every kernel is a Warp-CUDA operation; no authoritative host state participates.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

import warp as wp

from aleph.components.motor.hand import bell_off_rate

__all__ = [
    "ERMSpectrumModel",
    "ERMSpectrumParams",
    "erm_weld_force_kernel",
    "erm_clutch_force_kernel",
    "erm_bond_density_commit_kernel",
    "OFF_RATE_BELL_SLIP",
    "OFF_RATE_PEREVERZEV_CATCH_SLIP",
]

# device selector for the model-3 off-rate law (no biological default is implied by either integer)
OFF_RATE_BELL_SLIP = 0
OFF_RATE_PEREVERZEV_CATCH_SLIP = 1


class ERMSpectrumModel(Enum):
    """The three ERM membrane--cortex coupling models, ordered rigid -> dynamic."""

    RIGID_WELD = "rigid_weld"
    ELASTIC_FRICTION = "elastic_friction"
    DYNAMIC_CLUTCH = "dynamic_clutch"


@dataclass(frozen=True, slots=True)
class ERMSpectrumParams:
    """Complete source-gated contract for one ERM connector population under a chosen spectrum model.

    Only the fields the selected ``model`` uses are required; each required field must be finite and positive
    and the whole contract must carry an explicit provenance ``source`` (Fehon-class review is motivation, not
    a single-molecule datum — the operating constants are still PI-gated).  No field has a biological default.

    Fields by model:
        RIGID_WELD:       ``k_weld_pn_um``, ``rest_length_um``.
        ELASTIC_FRICTION: ``k_b_pn_um``, ``rho_b``, ``xi_pn_s_um``, ``rest_length_um``.
        DYNAMIC_CLUTCH:   ELASTIC_FRICTION fields + ``k_on_s``, ``rho_max``, and the off-rate law
                          (``off_rate_mode``): Bell -> ``k_off0_s`` + ``bell_force_pn``; Pereverzev ->
                          ``k_catch0_s`` + ``x_catch_um`` + ``k_slip0_s`` + ``x_slip_um`` + ``kT_pn_um``.
    """

    model: ERMSpectrumModel
    source: str
    # geometry / mechanics
    rest_length_um: float = 0.0
    k_weld_pn_um: float | None = None
    k_b_pn_um: float | None = None
    rho_b: float | None = None
    xi_pn_s_um: float | None = None
    # dynamic-clutch binding kinetics
    k_on_s: float | None = None
    rho_max: float | None = None
    off_rate_mode: int = OFF_RATE_BELL_SLIP
    k_off0_s: float | None = None
    bell_force_pn: float | None = None
    k_catch0_s: float | None = None
    x_catch_um: float | None = None
    k_slip0_s: float | None = None
    x_slip_um: float | None = None
    kT_pn_um: float | None = None

    def __post_init__(self) -> None:
        if not self.source.strip():
            raise ValueError("ERM spectrum params require an explicit literature/provenance source")
        if not math.isfinite(self.rest_length_um) or self.rest_length_um < 0.0:
            raise ValueError("rest_length_um must be finite and nonnegative")
        if self.model is ERMSpectrumModel.DYNAMIC_CLUTCH:
            if self.off_rate_mode not in (OFF_RATE_BELL_SLIP, OFF_RATE_PEREVERZEV_CATCH_SLIP):
                raise ValueError("off_rate_mode must be OFF_RATE_BELL_SLIP or OFF_RATE_PEREVERZEV_CATCH_SLIP")
        required = self._required_positive_fields()
        for name in required:
            value = getattr(self, name)
            if value is None or not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"model {self.model.value} requires finite positive {name}")
        if self.model is ERMSpectrumModel.DYNAMIC_CLUTCH:
            if self.rho_max is not None and self.rho_b is not None and self.rho_b > self.rho_max:
                raise ValueError("initial rho_b must not exceed rho_max")

    def _required_positive_fields(self) -> tuple[str, ...]:
        if self.model is ERMSpectrumModel.RIGID_WELD:
            return ("k_weld_pn_um",)
        base = ("k_b_pn_um", "rho_b", "xi_pn_s_um")
        if self.model is ERMSpectrumModel.ELASTIC_FRICTION:
            return base
        # DYNAMIC_CLUTCH
        kinetics = ("k_on_s", "rho_max")
        if self.off_rate_mode == OFF_RATE_BELL_SLIP:
            law = ("k_off0_s", "bell_force_pn")
        else:
            law = ("k_catch0_s", "x_catch_um", "k_slip0_s", "x_slip_um", "kT_pn_um")
        return base + kinetics + law


@wp.func
def _erm_catch_slip_off_rate(
    f: wp.float64, kc0: wp.float64, xc: wp.float64, ks0: wp.float64, xs: wp.float64, kT: wp.float64,
) -> wp.float64:
    """Pereverzev two-pathway catch/slip off-rate ``kc0·exp(-f·xc/kT) + ks0·exp(+f·xs/kT)`` [1/s].

    Matches the host oracle ``erm_spectrum_analytic.erm_catch_slip_off_rate`` and the FF clutch off-rate
    form; ``f`` is the per-bond load magnitude.
    """
    fa = wp.abs(f)
    return kc0 * wp.exp(-fa * xc / kT) + ks0 * wp.exp(fa * xs / kT)


@wp.kernel
def erm_weld_force_kernel(
    pos: wp.array(dtype=wp.vec3d),
    membrane_idx: wp.array(dtype=wp.int32),
    cortex_idx: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    k_weld: wp.float64,
    rest: wp.array(dtype=wp.float64),
    force: wp.array(dtype=wp.vec3d),
) -> None:
    """Model 1 (완전결합): accumulate the bilateral rigid-weld penalty ``k_weld·(L - rest)``.

    Bilateral — extension pulls the pair together, compression pushes it apart — so the connector enforces the
    Δu = 0 co-motion constraint as ``k_weld`` grows.  This is NOT the physical unilateral tether; it is the
    explicit infinitely-strong-attachment reference limit.  State is never mutated.
    """
    k = wp.tid()
    if bound[k] == 0:
        return
    mi = membrane_idx[k]
    ci = cortex_idx[k]
    delta = pos[mi] - pos[ci]
    length = wp.length(delta)
    if length < wp.float64(1.0e-12):
        return
    tension = k_weld * (length - rest[k])
    weld_force = tension * delta / length
    wp.atomic_add(force, mi, -weld_force)
    wp.atomic_add(force, ci, weld_force)


@wp.kernel
def erm_clutch_force_kernel(
    pos: wp.array(dtype=wp.vec3d),
    vel: wp.array(dtype=wp.vec3d),
    membrane_idx: wp.array(dtype=wp.int32),
    cortex_idx: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    k_b: wp.float64,
    rho_b: wp.array(dtype=wp.float64),
    xi: wp.float64,
    rest: wp.array(dtype=wp.float64),
    force: wp.array(dtype=wp.vec3d),
) -> None:
    """Models 2 & 3 (탄성·마찰 / 동적결합): unilateral density-scaled tension + along-bond viscous friction.

    Spring branch ``rho_b·k_b·max(L-rest,0)`` is unilateral (ezrin bears tension, slides, never pushes) and
    scaled by the local bond occupancy ``rho_b[k]`` — a fixed array for model 2, the state advanced by
    :func:`erm_bond_density_commit_kernel` for model 3.  Friction branch ``xi·du_dot`` acts along the bond
    (``du_dot`` = ``(v_m - v_c)·n̂``), opposing relative sliding in either direction; an overdamped/static
    caller passes a zero ``vel`` array to disable it.  State is never mutated.
    """
    k = wp.tid()
    if bound[k] == 0:
        return
    mi = membrane_idx[k]
    ci = cortex_idx[k]
    delta = pos[mi] - pos[ci]
    length = wp.length(delta)
    if length < wp.float64(1.0e-12):
        return
    n_hat = delta / length
    ext = length - rest[k]
    spring = wp.float64(0.0)
    if ext > wp.float64(0.0):
        spring = rho_b[k] * k_b * ext
    du_dot = wp.dot(vel[mi] - vel[ci], n_hat)
    axial = spring + xi * du_dot
    clutch_force = axial * n_hat
    wp.atomic_add(force, mi, -clutch_force)
    wp.atomic_add(force, ci, clutch_force)


@wp.kernel
def erm_bond_density_commit_kernel(
    pos: wp.array(dtype=wp.vec3d),
    membrane_idx: wp.array(dtype=wp.int32),
    cortex_idx: wp.array(dtype=wp.int32),
    bound: wp.array(dtype=wp.int32),
    k_b: wp.float64,
    rest: wp.array(dtype=wp.float64),
    rho_b: wp.array(dtype=wp.float64),
    k_on: wp.float64,
    rho_max: wp.float64,
    off_rate_mode: wp.int32,
    k_off0: wp.float64,
    bell_force: wp.float64,
    kc0: wp.float64,
    xc: wp.float64,
    ks0: wp.float64,
    xs: wp.float64,
    kT: wp.float64,
    tau: wp.float64,
    accepted: wp.array(dtype=wp.int32),
) -> None:
    """Model 3 (동적결합): advance the mean-field bond density at an accepted outer physical-time boundary.

    Per-bond load ``f = k_b·max(L-rest,0)`` sets the off-rate (Bell slip or Pereverzev catch/slip per
    ``off_rate_mode``); the linear balance ``d(rho)/dt = k_on·(rho_max - rho) - k_off·rho`` is advanced by the
    EXACT constant-coefficient tick ``rho <- rho_ss + (rho - rho_ss)·exp(-(k_on+k_off)·tau)`` — unconditionally
    stable and monotone toward ``rho_ss = k_on·rho_max/(k_on+k_off)`` for any ``tau``.  The device ``accepted``
    latch makes a rejected transaction a bit-exact no-op with no host decision.  Force is never touched here.
    """
    k = wp.tid()
    if accepted[0] == 0 or bound[k] == 0:
        return
    length = wp.length(pos[membrane_idx[k]] - pos[cortex_idx[k]])
    ext = wp.max(length - rest[k], wp.float64(0.0))
    f = k_b * ext
    if off_rate_mode == OFF_RATE_PEREVERZEV_CATCH_SLIP:
        k_off = _erm_catch_slip_off_rate(f, kc0, xc, ks0, xs, kT)
    else:
        k_off = bell_off_rate(f, k_off0, bell_force)
    rho_ss = k_on * rho_max / (k_on + k_off)
    decay = wp.exp(-(k_on + k_off) * tau)
    rho_new = rho_ss + (rho_b[k] - rho_ss) * decay
    rho_new = wp.max(wp.float64(0.0), wp.min(rho_new, rho_max))
    rho_b[k] = rho_new
