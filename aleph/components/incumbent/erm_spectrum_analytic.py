"""Host acceptance oracles for the 3-model ERM membrane--cortex connector spectrum.

The runtime CUDA primitives live in :mod:`aleph.components.incumbent.erm_spectrum`.  This module is acceptance
algebra only: it defines each model's closed-form axial force and (for the dynamic clutch) the bond-density
kinetics without advancing any biological state.  It is the ground truth the device kernels reproduce and the
only part of the spectrum that is exercised on the CUDA-free dev Mac.

Motivation (Fehon, McClatchey & Bretscher 2010, *Nat Rev Mol Cell Biol* 11:276, doi:10.1038/nrm2866): ERM
proteins are NOT a static membrane--cortex weld.  They are conformationally auto-inhibited linkers that are
switched on (PIP2 + phospho-activation), bind F-actin, bear load, and turn over — a **regulated dynamic clutch**
whose local density sets how strongly the plasma membrane is pinned to the cortex.  A faithful connector
therefore spans a spectrum whose limits bracket that biology:

* **model 1 -- RIGID_WELD (완전결합).**  The Δu = 0 co-motion constraint, realized as a *bilateral* stiff
  penalty ``F = k_weld·(L - r)`` (restores both extension and compression).  As ``k_weld -> inf`` the membrane
  node co-moves with its cortex node.  This is the infinitely-strong-attachment limit — not physical ERM, but
  the reference a clutch must reduce to when fully engaged and stiff.
* **model 2 -- ELASTIC_FRICTION (탄성·마찰).**  A Kelvin--Voigt clutch of *fixed* bond occupancy:
  ``F = rho_b·k_b·max(L - r, 0) + xi·du_dot``.  The spring is **unilateral** — ezrin bears tension and slides
  but cannot push the surfaces apart (Korkmazhan & Dunn 2022, Sci Adv 8:eabm1174) — scaled by the standing bond
  density ``rho_b``; ``xi·du_dot`` is the viscous membrane--cortex slip friction along the bond.  ``rho_b`` is a
  parameter here, not a state.
* **model 3 -- DYNAMIC_CLUTCH (동적결합).**  Model 2's force with an *evolving* occupancy obeying the mean-field
  clutch balance ``d(rho_b)/dt = k_on·(rho_max - rho_b) - k_off(f)·rho_b`` on the per-bond load
  ``f = k_b·max(L - r, 0)``, with ``k_off(f)`` either Bell slip or Pereverzev catch/slip.  This is the
  continuum mean-field generalization of the discrete per-tether Bell KMC in :mod:`erm_tether` (the
  ``rho_max = 1``, single-linker limit): membrane--cortex coupling strength emerges from binding vs
  load-dependent unbinding, exactly Fehon's dynamic-clutch homeostat.

Sanity Gate:
    * Dimensions: ``k [pN/um] * extension [um] = tension [pN]``; ``xi [pN·s/um] * du_dot [um/s] = [pN]``;
      ``k_on, k_off [1/s]``; ``rho_b, rho_max`` dimensionless bond occupancies.
    * Boundary cases: model-2/3 unilateral springs give exactly zero spring tension in compression and at rest;
      model-1 is bilateral (nonzero, sign-flipped, under compression).
    * Sign sense: extension cannot lower a tension-only spring's load; catch/slip off-rate falls then rises.
    * Kinetics limits: rho -> rho_max as k_off -> 0; rho -> 0 as k_off -> inf; the exact constant-coefficient
      tick relaxes toward the steady state monotonically.
    * Spectrum continuity: model 3 with ``k_off = 0`` and rho pinned at ``rho_max`` reproduces model 2's
      spring; model 2 with ``rho_b·k_b = k_weld`` and no friction reproduces model 1's tensile branch.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from aleph.laws.hand_kmc import pereverzev_off_rate

__all__ = [
    "erm_rigid_weld_axial_force",
    "erm_elastic_friction_axial_force",
    "erm_bell_slip_off_rate",
    "erm_catch_slip_off_rate",
    "erm_catch_slip_peak_force",
    "erm_bond_density_rhs",
    "erm_bond_density_steady_state",
    "erm_bond_density_tick",
]


def _validate_stiffness(name: str, value: float) -> None:
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")


def _as_length(length: npt.ArrayLike) -> npt.NDArray[np.float64]:
    values = np.asarray(length, dtype=np.float64)
    if np.any(~np.isfinite(values)) or np.any(values < 0.0):
        raise ValueError("length must contain finite nonnegative values")
    return values


def erm_rigid_weld_axial_force(
    length: npt.ArrayLike,
    *,
    k_weld: float,
    rest_length: float,
) -> npt.NDArray[np.float64]:
    """Model 1 bilateral weld axial force ``k_weld·(L - rest)`` in pN (signed: + = tensile).

    Bilateral because a rigid co-motion constraint resists both extension (pull together) and compression
    (push apart); as ``k_weld -> inf`` the relative displacement ``L - rest -> 0``.  Contrast the unilateral
    tension-only ERM spring of models 2/3.
    """
    _validate_stiffness("k_weld", k_weld)
    if not np.isfinite(rest_length) or rest_length < 0.0:
        raise ValueError("rest_length must be finite and nonnegative")
    values = _as_length(length)
    return np.asarray(k_weld * (values - rest_length), dtype=np.float64)


def erm_elastic_friction_axial_force(
    length: npt.ArrayLike,
    du_dot: npt.ArrayLike,
    *,
    k_b: float,
    rho_b: npt.ArrayLike,
    xi: float,
    rest_length: float,
) -> npt.NDArray[np.float64]:
    """Model 2 Kelvin--Voigt clutch axial force ``rho_b·k_b·max(L-rest,0) + xi·du_dot`` in pN.

    The spring branch is unilateral (tension only, ezrin cannot push) and scaled by the bond occupancy
    ``rho_b``; the friction branch ``xi·du_dot`` is a viscous along-bond dashpot that opposes relative sliding
    in either direction (``du_dot`` = along-bond relative speed ``(v_m - v_c)·n̂`` in um/s).
    """
    _validate_stiffness("k_b", k_b)
    if not np.isfinite(xi) or xi < 0.0:
        raise ValueError("xi must be finite and nonnegative")
    if not np.isfinite(rest_length) or rest_length < 0.0:
        raise ValueError("rest_length must be finite and nonnegative")
    values = _as_length(length)
    rate = np.asarray(du_dot, dtype=np.float64)
    occupancy = np.asarray(rho_b, dtype=np.float64)
    if np.any(~np.isfinite(rate)):
        raise ValueError("du_dot must be finite")
    if np.any(~np.isfinite(occupancy)) or np.any(occupancy < 0.0):
        raise ValueError("rho_b must be finite and nonnegative")
    spring = occupancy * k_b * np.maximum(values - rest_length, 0.0)
    return np.asarray(spring + xi * rate, dtype=np.float64)


def erm_bell_slip_off_rate(
    f: npt.ArrayLike,
    *,
    k_off0: float,
    bell_force: float,
) -> npt.NDArray[np.float64]:
    """Bell 1978 slip off-rate ``k_off0·exp(|f|/bell_force)`` in 1/s (mirrors the device ``bell_off_rate``)."""
    _validate_stiffness("k_off0", k_off0)
    _validate_stiffness("bell_force", bell_force)
    load = np.asarray(f, dtype=np.float64)
    if np.any(~np.isfinite(load)):
        raise ValueError("f must be finite")
    return np.asarray(k_off0 * np.exp(np.abs(load) / bell_force), dtype=np.float64)


def erm_catch_slip_off_rate(
    f: npt.ArrayLike,
    *,
    k_catch0: float,
    x_catch_um: float,
    k_slip0: float,
    x_slip_um: float,
    kT: float,
) -> npt.NDArray[np.float64]:
    """Pereverzev two-pathway catch/slip off-rate in 1/s (delegates to the FF single source of truth).

    ``k_catch0·exp(-f·x_catch/kT) + k_slip0·exp(+f·x_slip/kT)``: falls with load (catch) up to a peak F*,
    then rises (slip).  Faithful ERM: single-molecule ezrin--F-actin bonds are catch bonds under moderate load.
    """
    for name, value in (("k_catch0", k_catch0), ("x_catch_um", x_catch_um),
                        ("k_slip0", k_slip0), ("x_slip_um", x_slip_um), ("kT", kT)):
        _validate_stiffness(name, value)
    load = np.asarray(f, dtype=np.float64)
    if np.any(~np.isfinite(load)):
        raise ValueError("f must be finite")
    return np.asarray(pereverzev_off_rate(load, k_catch0, x_catch_um, k_slip0, x_slip_um, kT),
                      dtype=np.float64)


def erm_catch_slip_peak_force(
    *,
    k_catch0: float,
    x_catch_um: float,
    k_slip0: float,
    x_slip_um: float,
    kT: float,
) -> float:
    """Analytic catch/slip peak load ``F* = kT/(x_c+x_s)·ln[(k_c0 x_c)/(k_s0 x_s)]`` in pN (min off-rate).

    Solving ``d(off_rate)/df = 0``.  Requires the catch branch to dominate at low load
    (``k_c0·x_c > k_s0·x_s``); otherwise the bond is slip-only and there is no catch peak.
    """
    for name, value in (("k_catch0", k_catch0), ("x_catch_um", x_catch_um),
                        ("k_slip0", k_slip0), ("x_slip_um", x_slip_um), ("kT", kT)):
        _validate_stiffness(name, value)
    ratio = (k_catch0 * x_catch_um) / (k_slip0 * x_slip_um)
    if ratio <= 1.0:
        raise ValueError("no catch peak: k_catch0·x_catch must exceed k_slip0·x_slip (slip-only bond)")
    return float(kT / (x_catch_um + x_slip_um) * np.log(ratio))


def erm_bond_density_rhs(
    rho_b: npt.ArrayLike,
    *,
    k_on: float,
    rho_max: float,
    k_off: npt.ArrayLike,
) -> npt.NDArray[np.float64]:
    """Model 3 mean-field clutch balance ``d(rho_b)/dt = k_on·(rho_max - rho_b) - k_off·rho_b`` in 1/s.

    ``k_off`` is the load-dependent off-rate evaluated at the per-bond load (Bell or Pereverzev); the caller
    supplies it so the same balance serves both catch/slip variants.
    """
    _validate_stiffness("k_on", k_on)
    _validate_stiffness("rho_max", rho_max)
    rho = np.asarray(rho_b, dtype=np.float64)
    off = np.asarray(k_off, dtype=np.float64)
    if np.any(~np.isfinite(rho)) or np.any(rho < 0.0):
        raise ValueError("rho_b must be finite and nonnegative")
    if np.any(~np.isfinite(off)) or np.any(off < 0.0):
        raise ValueError("k_off must be finite and nonnegative")
    return np.asarray(k_on * (rho_max - rho) - off * rho, dtype=np.float64)


def erm_bond_density_steady_state(
    *,
    k_on: float,
    rho_max: float,
    k_off: npt.ArrayLike,
) -> npt.NDArray[np.float64]:
    """Steady-state occupancy ``rho_ss = k_on·rho_max/(k_on + k_off)`` (root of the model-3 balance)."""
    _validate_stiffness("k_on", k_on)
    _validate_stiffness("rho_max", rho_max)
    off = np.asarray(k_off, dtype=np.float64)
    if np.any(~np.isfinite(off)) or np.any(off < 0.0):
        raise ValueError("k_off must be finite and nonnegative")
    return np.asarray(k_on * rho_max / (k_on + off), dtype=np.float64)


def erm_bond_density_tick(
    rho_b: npt.ArrayLike,
    *,
    k_on: float,
    rho_max: float,
    k_off: npt.ArrayLike,
    tau: float,
) -> npt.NDArray[np.float64]:
    """Exact constant-coefficient advance of the model-3 balance over one accepted tick ``tau`` [s].

    Freezing ``k_off`` over the tick, the linear ODE integrates exactly to
    ``rho(t+tau) = rho_ss + (rho - rho_ss)·exp(-(k_on + k_off)·tau)`` — unconditionally stable and monotone
    toward ``rho_ss`` for any ``tau``, unlike an explicit Euler step which can overshoot.  This is the host
    reference for the device commit kernel; the kernel advances only on an accepted outer step.
    """
    if not np.isfinite(tau) or tau < 0.0:
        raise ValueError("tau must be finite and nonnegative")
    rho = np.asarray(rho_b, dtype=np.float64)
    if np.any(~np.isfinite(rho)) or np.any(rho < 0.0):
        raise ValueError("rho_b must be finite and nonnegative")
    off = np.asarray(k_off, dtype=np.float64)
    if np.any(~np.isfinite(off)) or np.any(off < 0.0):
        raise ValueError("k_off must be finite and nonnegative")
    rho_ss = erm_bond_density_steady_state(k_on=k_on, rho_max=rho_max, k_off=off)
    decay = np.exp(-(k_on + off) * tau)
    return np.asarray(rho_ss + (rho - rho_ss) * decay, dtype=np.float64)
