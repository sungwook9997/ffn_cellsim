r"""Bell slip + Pereverzev catch-slip detachment kinetics — analytic ground truth for the per-head off-rate.

Host-side acceptance oracle (pure NumPy, NO Warp): the Warp-CUDA per-head detachment kernels and the device
engaged fraction are gated AGAINST these closed forms.

TWO force-dependent off-rate laws live here, one per bond CLASS — both mechanistic Bell-Evans, NOT a
detailed-balance Metropolis proxy (CLAUDE.md worked example: "Bond off-rate: Bell-Evans force-dependent"):

  1. Pure Bell SLIP (:func:`bell_off_rate`) — off-rate RISES monotonically with load:

        p_off(f) = k_off0 * exp(|f| / f0),   f0 = kBT / x_beta               (Bell 1978 slip)

     This is the correct law for the passive, non-processive linkers: alpha-actinin, ERM/ezrin, kinesin,
     and the reforming crosslinkers (NF2007 §10.1). ``k_off0`` is the zero-force rate, ``f0`` the Bell
     characteristic force set by the bond length ``x_beta``.

  2. Pereverzev two-pathway CATCH-SLIP (:func:`catch_slip_off_rate`) — the off-rate FALLS with load up to a
     peak-lifetime force ``F*`` (catch), then RISES (slip):

        p_off(f) = k_catch0 * exp(-|f| * x_catch / kBT) + k_slip0 * exp(+|f| * x_slip / kBT)

     ⚠ FIDELITY CORRECTION (2026-07-23, sourcing §b): the NMII HEAD-ACTIN bond is catch-slip, NOT pure
     slip. Kovacs, Thirumurugan, Knight & Sellers 2007, *PNAS* 104(24):9994 (10.1073/pnas.0701181104)
     show resistive load RAISES the bound-head duty for tension maintenance — ADP release (the
     rate-limiting detachment step) is slowed 5x (NM2A) / 12x (NM2B) under load. A pure Bell slip has the
     OPPOSITE sign (duty falls under load); modelling the NMII head as slip understates sustained cortical
     tension. The catch pathway (``k_catch0``, ``x_catch``) is that load-slowed ADP-release detachment; the
     slip pathway (``k_slip0``, ``x_slip``) is forced unbinding at high load. Same two-pathway form the
     project already uses for filamin/integrin/ERM (``ff.hand_kmc.pereverzev_off_rate``, the single source
     of truth this module delegates to). The pure-slip law above is UNCHANGED — it is still correct for the
     genuine slip bonds and is what alpha-actinin/ERM/crosslink/kinesin consumers import.

The engaged fraction is EMERGENT, not imposed (P2 "Bell-emergent engaged fraction"): a head cycles free
<-> bound with attach rate ``k_on`` and detach rate ``p_off(f)``, so its steady-state bound probability is

    phi_b(f) = k_on / (k_on + p_off(f))                                       (two-state steady state)

For a SLIP bond (:func:`engaged_fraction_steady`) this FALLS with load (heads shed). For the NMII CATCH-SLIP
bond (:func:`engaged_fraction_catch_slip`) it RISES to a maximum at ``F*`` then falls — the Kovacs 2007
load-borne tension-maintenance behaviour. At zero load ``phi_b(0) = k_on/(k_on + p_off(0))`` reconciles with
the measured unloaded duty (NM2B ~0.2-0.3, Nagy 2013 / Kovacs 2003) — a cross-check, never a runtime input.
The ensemble stall then emerges as (engaged count) x (per-head stall), NOT as ``N_side * F_head``.

Per-tick KMC probabilities (Poisson, ported from ff/hand_kmc.py so the device kernel matches the oracle):

    P_attach = 1 - exp(-tau * k_on),   P_detach = 1 - exp(-tau * p_off)

Sanity Gate (of the oracle itself — tests/ac/motor/test_bell_kinetics_oracle.py):
  * SLIP: p_off(0) = k_off0; p_off strictly increasing in |f|; symmetric in sign of f;
  * f0 = kBT/x_beta reproduces the ff/hand_kmc NMIIA value (~7.1 pN at x_beta=0.6 nm);
  * CATCH-SLIP: p_off falls then rises (biphasic), minimum at the analytic F*; reduces to pure slip as the
    catch prefactor -> 0; matches the ff.hand_kmc.pereverzev_off_rate SoT bit-for-formula;
  * probabilities in [0, 1); small-tau expansion P ~ tau*rate; P -> 1 as tau -> inf;
  * phi_b in (0, 1]; SLIP phi_b strictly decreasing in load; CATCH-SLIP phi_b RISES then falls, peak at F*,
    phi_b(F*) > phi_b(0) (load raises duty — the Kovacs 2007 sign);
  * detailed-balance sanity: with k_on -> inf, phi_b -> 1 (always bound); k_on -> 0, phi_b -> 0.

References:
  Bell, G.I. 1978, Science 200:618 (slip bond); Evans 1997. Pereverzev et al. 2005, Biophys J 89:1446
  (two-pathway catch-slip). Kovacs, Thirumurugan, Knight & Sellers 2007, PNAS 104(24):9994
  (10.1073/pnas.0701181104 — NMII duty RISES under load). Nagy 2013 JBC 10.1074/jbc.M112.424671 (NM2B duty
  0.2-0.3). Nedelec & Foethke 2007 §10.1 (Cytosim Hand). Veigel 2002 (x_beta).
"""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt

from aleph.laws.hand_kmc import pereverzev_off_rate as _pereverzev_off_rate

__all__ = [
    "bell_off_rate",
    "bell_f0_from_x_beta",
    "catch_slip_off_rate",
    "catch_slip_peak_force",
    "attach_probability",
    "detach_probability",
    "engaged_fraction_steady",
    "engaged_fraction_catch_slip",
    "bound_state_lifetime",
]

# thermal scale at 37 C (pN*um); matches ff/units KBT and params_i0b3.yaml::kBT.
KBT_PN_UM = 4.28e-3


def bell_off_rate(
    force: npt.ArrayLike,
    k_off0: float,
    f0: float,
) -> npt.NDArray[np.float64]:
    """Bell slip off-rate ``p_off = k_off0 * exp(|f| / f0)`` [1/s].

    Args:
        force: load magnitude ``f`` [pN] on the bound head (scalar or array; sign ignored).
        k_off0: zero-force detachment rate [1/s] (> 0).
        f0: Bell characteristic force [pN] (> 0), ``f0 = kBT / x_beta``.

    Returns:
        Force-dependent off-rate [1/s], same shape as ``force``.
    """
    if k_off0 <= 0.0:
        raise ValueError("k_off0 must be positive")
    if f0 <= 0.0:
        raise ValueError("f0 must be positive")
    f = np.asarray(force, dtype=np.float64)
    return k_off0 * np.exp(np.abs(f) / f0)


def bell_f0_from_x_beta(x_beta_um: float, kBT: float = KBT_PN_UM) -> float:
    """Bell characteristic force ``f0 = kBT / x_beta`` [pN] from bond length ``x_beta`` [µm]."""
    if x_beta_um <= 0.0:
        raise ValueError("x_beta must be positive")
    return kBT / x_beta_um


def catch_slip_off_rate(
    force: npt.ArrayLike,
    k_catch0: float,
    x_catch: float,
    k_slip0: float,
    x_slip: float,
    kT: float = KBT_PN_UM,
) -> npt.NDArray[np.float64]:
    r"""Pereverzev two-pathway catch-slip off-rate for the NMII head-actin bond [1/s].

    ``p_off(f) = k_catch0 * exp(-|f|*x_catch/kT) + k_slip0 * exp(+|f|*x_slip/kT)``.

    The CATCH pathway (``k_catch0``, ``x_catch``) is the load-slowed ADP-release detachment (Kovacs 2007):
    it makes the off-rate FALL with load up to the peak-lifetime force ``F*`` (see :func:`catch_slip_peak_force`),
    so the head stays bound LONGER under resistive load — the opposite sign of a pure Bell slip. The SLIP
    pathway (``k_slip0``, ``x_slip``) is forced unbinding that dominates past ``F*``, so the off-rate rises
    again (biphasic). With ``k_catch0 -> 0`` this reduces to a pure Bell slip ``k_slip0*exp(|f|*x_slip/kT)``.

    Load magnitude is used (``|f|``); the sign of ``force`` is ignored, matching :func:`bell_off_rate` and the
    device ``wp.abs(f)``. Delegates to :func:`aleph.laws.hand_kmc.pereverzev_off_rate` — the project single
    source of truth for this two-pathway form (filamin / integrin / ERM all reuse it).

    Args:
        force: load magnitude ``f`` [pN] on the bound head (scalar or array; sign ignored).
        k_catch0: zero-force catch-pathway rate [1/s] (> 0).
        x_catch: catch-pathway bond length [µm] (> 0); larger => stronger load-strengthening.
        k_slip0: zero-force slip-pathway rate [1/s] (> 0).
        x_slip: slip-pathway (Bell) bond length [µm] (> 0).
        kT: thermal energy [pN·µm] (default 37 °C engine value).

    Returns:
        Catch-slip off-rate [1/s], same shape as ``force``.
    """
    for name, value in (("k_catch0", k_catch0), ("x_catch", x_catch),
                        ("k_slip0", k_slip0), ("x_slip", x_slip), ("kT", kT)):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    f = np.abs(np.asarray(force, dtype=np.float64))
    return np.asarray(_pereverzev_off_rate(f, k_catch0, x_catch, k_slip0, x_slip, kT), dtype=np.float64)


def catch_slip_peak_force(
    k_catch0: float,
    x_catch: float,
    k_slip0: float,
    x_slip: float,
    kT: float = KBT_PN_UM,
) -> float:
    r"""Catch-slip peak-lifetime load ``F* = kT/(x_catch+x_slip) * ln[(k_catch0*x_catch)/(k_slip0*x_slip)]`` [pN].

    The load where ``d(p_off)/df = 0`` — the off-rate MINIMUM, hence the bound-lifetime and engaged-fraction
    MAXIMUM. This is the Kovacs-2007 tension-maintenance operating load: duty rises up to ``F*`` then falls.

    Requires the catch branch to dominate at low load (``k_catch0*x_catch > k_slip0*x_slip``); otherwise the
    bond is slip-only (no interior peak) and this raises ``ValueError``.
    """
    for name, value in (("k_catch0", k_catch0), ("x_catch", x_catch),
                        ("k_slip0", k_slip0), ("x_slip", x_slip), ("kT", kT)):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    ratio = (k_catch0 * x_catch) / (k_slip0 * x_slip)
    if ratio <= 1.0:
        raise ValueError(
            "no interior catch peak: k_catch0*x_catch <= k_slip0*x_slip (bond is slip-dominated everywhere)"
        )
    return float(kT / (x_catch + x_slip) * np.log(ratio))


def attach_probability(tau: float, k_on: float) -> float:
    """Per-tick attach probability ``1 - exp(-tau * k_on)`` (Poisson; ~ tau*k_on for small tau*k_on)."""
    if tau < 0.0 or k_on < 0.0:
        raise ValueError("tau and k_on must be non-negative")
    return float(-np.expm1(-tau * k_on))


def detach_probability(tau: float, p_off: npt.ArrayLike) -> npt.NDArray[np.float64]:
    """Per-tick detach probability ``1 - exp(-tau * p_off)`` for off-rate ``p_off`` [1/s]."""
    if tau < 0.0:
        raise ValueError("tau must be non-negative")
    return -np.expm1(-tau * np.asarray(p_off, dtype=np.float64))


def engaged_fraction_steady(
    k_on: float,
    force: npt.ArrayLike,
    k_off0: float,
    f0: float,
) -> npt.NDArray[np.float64]:
    """Steady-state bound probability ``phi_b(f) = k_on / (k_on + p_off(f))`` of one head under load ``force``.

    The EMERGENT engaged fraction (P2): it falls with load because the Bell slip off-rate rises. At zero
    load it is ``k_on / (k_on + k_off0)`` (the unloaded duty cross-check). This is the mean-field bound
    probability; the stochastic head count is Binomial(N_side, phi_b) — see ``ensemble_stall_analytic``.

    Args:
        k_on: per-head attach rate [1/s] (>= 0).
        force: per-head load [pN] (scalar or array).
        k_off0: zero-force off-rate [1/s].
        f0: Bell characteristic force [pN].

    Returns:
        Steady-state bound probability in (0, 1], same shape as ``force``.
    """
    if k_on < 0.0:
        raise ValueError("k_on must be non-negative")
    p_off = bell_off_rate(force, k_off0, f0)
    return k_on / (k_on + p_off)


def engaged_fraction_catch_slip(
    k_on: float,
    force: npt.ArrayLike,
    k_catch0: float,
    x_catch: float,
    k_slip0: float,
    x_slip: float,
    kT: float = KBT_PN_UM,
) -> npt.NDArray[np.float64]:
    r"""Steady-state bound probability ``phi_b(f) = k_on / (k_on + p_off(f))`` for the NMII CATCH-SLIP head.

    Unlike the pure-slip :func:`engaged_fraction_steady` (which falls monotonically with load), this RISES
    with load up to the peak-lifetime force ``F*`` (:func:`catch_slip_peak_force`), then falls — the Kovacs
    2007 sign, where resistive load raises the bound duty for tension maintenance. At zero load
    ``phi_b(0) = k_on / (k_on + k_catch0 + k_slip0)`` is the unloaded duty (NM2B ~0.2-0.3 cross-check).

    Args:
        k_on: per-head attach rate [1/s] (>= 0).
        force: per-head load [pN] (scalar or array; sign ignored).
        k_catch0, x_catch, k_slip0, x_slip, kT: catch-slip law parameters (see :func:`catch_slip_off_rate`).

    Returns:
        Steady-state bound probability in (0, 1], same shape as ``force``.
    """
    if k_on < 0.0:
        raise ValueError("k_on must be non-negative")
    p_off = catch_slip_off_rate(force, k_catch0, x_catch, k_slip0, x_slip, kT)
    return k_on / (k_on + p_off)


def bound_state_lifetime(
    loads_pn: npt.ArrayLike,
    *,
    catch_slip: bool,
    k_off0: float | None = None,
    f0: float | None = None,
    k_catch0: float | None = None,
    x_catch: float | None = None,
    k_slip0: float | None = None,
    x_slip: float | None = None,
    kT: float = KBT_PN_UM,
) -> dict[str, float | int | str]:
    r"""Return the bound-state lifetime [s] a run's OWN kinetics set, with its derivation.

    WHY THIS IS SHARED CODE.  Until 2026-07-28 the cortex-motor driver carried this as a literal,
    ``1/0.4``, and all three numbers in that expression were wrong for the runs it stamped: this
    repository's Bell prefactor is 0.35, the runs used the CATCH-SLIP path whose zero-load off-rate is
    the SUM of the two pathways and not the Bell prefactor at all, and the heads were under load, which
    for a catch bond makes them live LONGER than either.  That is the denominator every stationarity
    contract is expressed in, so getting it wrong mis-scales the required window.  It lives here — with
    the closed forms it evaluates — rather than in a driver, because a driver-local copy is a defect
    each new driver re-acquires; that is exactly how the literal survived as long as it did.

    TWO VALUES, and the caller should multiply the LARGER by its declared window count:

    * ``zero_load_s`` is computable from the declared parameters BEFORE the run, so a contract can name
      it without having seen data.
    * ``mean_load_s`` is ``1/<k_off>`` over the bonds actually bound.  Measured, and for a CATCH bond it
      is the physically relevant one because load slows detachment.

    Taking the SHORTER of the two would let a loaded run declare settlement sooner — gate-loosening by
    arithmetic — so ``conservative_s`` is always the larger finite value.

    Args:
        loads_pn: Per-bond load magnitudes [pN] over the BOUND bonds only.  Empty is allowed and yields
            a NaN ``mean_load_s`` rather than a zero, because "nothing bound" is not "zero lifetime".
        catch_slip: Whether the run uses the Pereverzev two-pathway form.  When ``False`` the pure Bell
            slip applies and ``k_off0``/``f0`` are required.
        k_off0: Bell slip prefactor [1/s] (required when ``catch_slip`` is ``False``).
        f0: Bell force scale [pN] (required when ``catch_slip`` is ``False``).
        k_catch0: Catch-pathway zero-force rate [1/s] (required when ``catch_slip``).
        x_catch: Catch-pathway bond length [µm] (required when ``catch_slip``).
        k_slip0: Slip-pathway zero-force rate [1/s] (required when ``catch_slip``).
        x_slip: Slip-pathway bond length [µm] (required when ``catch_slip``).
        kT: Thermal energy [pN·µm].

    Returns:
        ``zero_load_s`` · ``mean_load_s`` · ``conservative_s`` · ``n_bound`` · ``derivation``.

    Raises:
        ValueError: If a rate needed by the selected pathway is missing or non-positive — a silently
            defaulted off-rate is the failure this function exists to prevent.
    """
    if catch_slip:
        for name, value in (("k_catch0", k_catch0), ("x_catch", x_catch),
                            ("k_slip0", k_slip0), ("x_slip", x_slip)):
            if value is None:
                raise ValueError(f"catch_slip lifetime needs {name}; refusing to default it")
        k_zero = float(k_catch0) + float(k_slip0)                       # type: ignore[arg-type]
        derivation = ("Pereverzev two-pathway at zero load: k_catch0 + k_slip0 = "
                      f"{k_catch0:g} + {k_slip0:g} = {k_zero:g} /s")
    else:
        if k_off0 is None or f0 is None:
            raise ValueError("Bell slip lifetime needs k_off0 and f0; refusing to default them")
        k_zero = float(k_off0)
        derivation = f"Bell slip prefactor k_off0 = {k_zero:g} /s at zero load"
    if not (k_zero > 0.0 and math.isfinite(k_zero)):
        raise ValueError(f"zero-load off-rate must be positive-finite; got {k_zero!r}")

    loads = np.abs(np.asarray(loads_pn, dtype=np.float64).ravel())
    n_bound = int(loads.size)
    if n_bound:
        if catch_slip:
            rates = catch_slip_off_rate(loads, float(k_catch0), float(x_catch),    # type: ignore[arg-type]
                                        float(k_slip0), float(x_slip), kT)         # type: ignore[arg-type]
        else:
            rates = bell_off_rate(loads, float(k_off0), float(f0))                 # type: ignore[arg-type]
        mean_rate = float(np.mean(rates))
        mean_load_s = 1.0 / mean_rate if mean_rate > 0.0 else float("nan")
    else:
        mean_load_s = float("nan")

    zero_load_s = 1.0 / k_zero
    finite = [v for v in (zero_load_s, mean_load_s) if math.isfinite(v)]
    return {
        "zero_load_s": zero_load_s,
        "mean_load_s": mean_load_s,
        "conservative_s": max(finite) if finite else float("nan"),
        "n_bound": n_bound,
        "derivation": (derivation + "; mean_load_s = 1/<k_off> over the bound bonds' own loads through "
                       "the same closed form the device kernel implements"),
    }
