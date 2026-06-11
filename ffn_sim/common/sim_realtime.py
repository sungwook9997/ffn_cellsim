"""Map ACCELERATED DCM simulation time → a REAL-TIME (biological) equivalent.

WHY THIS EXISTS (PI 2026-06-11)
-------------------------------
Every production run and figure currently shows only the raw integrator step
count, which the PI cannot compare to real experiments (those run over
minutes/hours/days). The DCM stack uses *accelerated kinetics*: the BAOAB
timestep is ``dt ~ 1e-9 s`` so the raw integrator time ``t_sim = dt · steps``
of a typical run is only ~µs — that is the MECHANICAL clock, NOT the biological
time. The protrusion / spreading machinery is run in explicit *kinetic
fast-forward*: every spreading rate is multiplied by a single, explicitly
reported acceleration factor ``S`` (see ``scripts/h7_spreading_dynamics.py`` and
``cell/spreading_drive.py``), so the EFFECTIVE (real-time-equivalent) clock runs
``S×`` faster than the mechanical clock:

    t_real  =  S · t_sim                                                   (1)

This module computes, for any ``(steps, dt)``:

    * ``t_sim``        = dt · steps                         [s]  (raw integrator)
    * ``t_real``       = S · t_sim                          [s]  (real-time equiv.)
    * ``accel_factor`` = S                                  [dimensionless]

HONESTY / CAVEAT (read this before trusting ``t_real``)
-------------------------------------------------------
``t_real`` is a **MAPPING** from accelerated-sim time to a real-time equivalent
*via the documented acceleration factor* — it is **NOT a native real-time
integration**. The mechanism is untouched (the BAOAB integrator still moves
every particle, the force-velocity shape of the protrusion is preserved); only
the absolute timescale of the rate-limiting kinetics is compressed by ``S``. So
``t_real`` answers "if the spreading front advanced at its real biological
velocity instead of the accelerated one, how long would this much front advance
take?" — it is a velocity-/rate-matched equivalent, not a measured wall-time and
not a second integration. Always display ``t_sim`` AND ``t_real`` side by side,
both clearly labelled (sim vs real), with the ``accel`` factor stated — never
show ``t_real`` alone as if it were native.

ACCELERATION FACTOR — DERIVATION (rate-limiting process)
--------------------------------------------------------
The rate-limiting biological process of cell spreading is **lamellipodial
protrusion / the advance of the cell-spreading front**. The model advances the
spread front at an effective (band-matched, accelerated) velocity ``v_sim``; the
real lamellipodial protrusion / cell-spreading-front velocity is ``v_real``. The
acceleration factor is the ratio of the two velocities:

    S = accel = v_sim / v_real                                             (2)

and then ``t_real = t_sim · S`` advances the real front the same DISTANCE the
accelerated front advances in ``t_sim`` — i.e. distance is the invariant, time is
rescaled by the velocity ratio. This is the same ``S`` that
``scripts/h7_spreading_dynamics.py`` reports as ``S_accel_effective`` (it backs
``S`` out of the emergent front advance ``Δr`` and the literature velocity:
``t_phys = Δr / v_real``, ``S = t_phys / t_sim``).

PRIMARY BASIS — spreading-front velocity (used by default)
..........................................................
* ``v_sim``  : the accelerated front velocity the protrusion engine drives. The
  H.7 dynamic-spreading run reached its in-band footprint (A/A₀ 1.0→2.03 ∈ [2,4])
  at an *effective* acceleration ``S ≈ 6e5`` (``reference_ffn_spreading_dynamics``
  / commit f654531): the ~µs run mapped to ~2.66 physiological MINUTES. With the
  literature mapping velocity ``v_real ≈ 1 µm/min`` (Betorz 2023 P1, ~3 µm in
  ~3 min) this fixes ``v_sim = S · v_real ≈ 6e5 µm/min`` for the DCM dt-band
  (dt=1e-9 s). We therefore anchor the default acceleration at ``S = 6e5``.
* ``v_real`` : real single-cell lamellipodial protrusion / spreading-front
  velocity. Literature band for the *leading-edge protrusion* rate is
  **~0.1–0.5 µm/min** (slow, persistent lamellipodial advance; Abercrombie 1970;
  Ponti 2004 actin retrograde-flow / protrusion; Giannone 2004 edge-advance
  cycles). The faster **whole-cell early-spreading front** band is **~3–12 µm/min**
  (Cuvelier 2007; Dubin-Thaler 2004; Betorz 2023). The default mapping uses the
  Betorz whole-cell value ``v_real = 1 µm/min`` so that the reported ``S`` matches
  the H.7 ``S_accel_effective ≈ 6e5`` exactly; the slow-protrusion band is exposed
  for the conservative (longer real-time) reading.

ALTERNATIVE BASIS — cell-cycle / proliferation rate
...................................................
For runs where DIVISION (not spreading) dominates the timescale, the rate-limiting
process is the **cell cycle**: real ~**12–24 h per division** (MCF7 doubling
~24–30 h; generic mammalian cycle 12–24 h). Here the invariant is *one division
event*: ``t_real_per_div = T_cycle`` and the acceleration is
``S_cycle = T_cycle / t_sim_per_div`` where ``t_sim_per_div = dt · div_every``.
``proliferation_accel`` exposes this mapping for division-dominated runs.

SANITY (against real MCF7 spreading)
------------------------------------
A typical DCM run of 12000 steps × 1e-9 s = 12 µs sim, at S=6e5, maps to
``t_real = 12e-6 · 6e5 = 7.2 s``… per *unit* — but spreading runs are reported by
the FRONT ADVANCE, and the H.7 anchor (the full in-band spread) lands at ~2.66
min, i.e. the minutes-scale of real early MCF7 spreading (which reaches a plateau
over ~minutes-to-tens-of-minutes, full spread over ~hours). The mapping puts a
~µs accelerated run into the seconds-to-minutes real band, scaling to hours for
the longer production runs — consistent with real MCF7 spreading taking minutes
(early) to hours (full). See ``RealTimeMapping`` for the per-run numbers.

UNITS: SI throughout (seconds, metres). Velocities accepted in µm/min for
convenience and converted internally.
"""

from __future__ import annotations

from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Literature constants (the rate-limiting-process bands; documented above)
# ---------------------------------------------------------------------------
#: Default kinetic acceleration factor S for the DCM spreading regime (dt~1e-9 s).
#: Anchored to the H.7 dynamic-spreading run's S_accel_effective ≈ 6e5
#: (reference_ffn_spreading_dynamics, commit f654531) — the run that reached the
#: in-band footprint A/A₀ 1.0→2.03 ∈ [2,4] and mapped to ~2.66 physiological min.
DEFAULT_SPREADING_ACCEL: float = 6.0e5

#: Real lamellipodial protrusion / cell-spreading-front velocity bands [µm/min].
#: Slow leading-edge protrusion (Abercrombie 1970; Ponti 2004; Giannone 2004):
V_REAL_PROTRUSION_UM_MIN_LO: float = 0.1
V_REAL_PROTRUSION_UM_MIN_HI: float = 0.5
#: Whole-cell early-spreading front (Cuvelier 2007; Dubin-Thaler 2004; Betorz 2023):
V_REAL_SPREAD_FRONT_UM_MIN_LO: float = 3.0
V_REAL_SPREAD_FRONT_UM_MIN_HI: float = 12.0
#: Betorz 2023 P1 single-cell mapping velocity used by h7_spreading_dynamics
#: (~3 µm advance in ~3 min) — the value that fixes S ≈ 6e5 at dt=1e-9 s.
V_REAL_BETORZ_UM_MIN: float = 1.0

#: Cell-cycle / division time band [s] (12–24 h; MCF7 doubling ~24–30 h).
T_CYCLE_S_LO: float = 12.0 * 3600.0
T_CYCLE_S_HI: float = 24.0 * 3600.0
T_CYCLE_S_DEFAULT: float = 24.0 * 3600.0  # MCF7-ish ~1 day per division

_UM_PER_MIN_TO_M_PER_S = 1.0e-6 / 60.0


@dataclass(frozen=True, slots=True)
class RealTimeMapping:
    """A computed sim→real time mapping for a single run (pure data).

    Attributes:
        steps: integrator steps run.
        dt_s: BAOAB timestep [s].
        t_sim_s: raw integrator time = dt · steps [s] (the mechanical clock, ~µs).
        accel_factor: the acceleration factor S used (dimensionless).
        t_real_s: real-time equivalent = S · t_sim [s].
        basis: which rate-limiting process the factor is derived from
            ("spreading_front" or "cell_cycle").
        detail: a one-line human description of the derivation (for figure
            captions / metrics provenance).
    """

    steps: int
    dt_s: float
    t_sim_s: float
    accel_factor: float
    t_real_s: float
    basis: str
    detail: str

    @property
    def t_sim_human(self) -> str:
        """Raw integrator time as a human-readable string (e.g. ``'12 µs'``)."""
        return format_time(self.t_sim_s)

    @property
    def t_real_human(self) -> str:
        """Real-time equivalent as a human-readable string (e.g. ``'7.2 s'``)."""
        return format_time(self.t_real_s)

    def to_dict(self) -> dict:
        """JSON-serialisable metrics block (drop straight into a metrics dict)."""
        return {
            "steps": int(self.steps),
            "dt_s": float(self.dt_s),
            "t_sim_s": float(self.t_sim_s),
            "accel_factor": float(self.accel_factor),
            "t_real_s": float(self.t_real_s),
            "t_sim_human": format_time(self.t_sim_s),
            "t_real_human": format_time(self.t_real_s),
            "realtime_basis": self.basis,
            "realtime_detail": self.detail,
        }

    def title_str(self) -> str:
        """Compact one-line readout for a figure/MP4 title.

        e.g. ``"t = 12 µs sim  ≈ 7.2 s real  (accel 6×10⁵, spreading_front)"``.
        """
        return (f"t = {format_time(self.t_sim_s)} sim  "
                f"≈ {format_time(self.t_real_s)} real  "
                f"(accel {_sci(self.accel_factor)}, {self.basis})")


def spreading_accel(*, v_sim_um_min: float | None = None,
                    v_real_um_min: float = V_REAL_BETORZ_UM_MIN,
                    accel: float | None = None) -> float:
    """Acceleration factor S from the spreading-front velocity ratio (eq. 2).

    ``S = v_sim / v_real``. Either pass an explicit ``accel`` (e.g. the run's
    reported ``S_accel_effective``), or pass the accelerated front velocity
    ``v_sim_um_min`` and the real velocity ``v_real_um_min`` and let the ratio be
    formed. With neither, returns ``DEFAULT_SPREADING_ACCEL`` (the H.7 anchor).

    Args:
        v_sim_um_min: accelerated (band-matched) front velocity [µm/min].
        v_real_um_min: real lamellipodial / spreading-front velocity [µm/min]
            (default = Betorz 2023 mapping value, 1 µm/min).
        accel: an explicit acceleration factor (overrides the velocity ratio).

    Returns:
        The dimensionless acceleration factor S (> 0).
    """
    if accel is not None:
        if accel <= 0:
            raise ValueError(f"accel must be > 0, got {accel}")
        return float(accel)
    if v_sim_um_min is not None:
        if v_real_um_min <= 0:
            raise ValueError(f"v_real_um_min must be > 0, got {v_real_um_min}")
        return float(v_sim_um_min) / float(v_real_um_min)
    return DEFAULT_SPREADING_ACCEL


def proliferation_accel(*, dt_s: float, div_every_steps: int,
                        t_cycle_s: float = T_CYCLE_S_DEFAULT) -> float:
    """Acceleration factor S for a DIVISION-dominated run (cell-cycle basis).

    The invariant is one division event: one accelerated division occurs every
    ``div_every_steps`` integration steps (= ``dt · div_every_steps`` of sim
    time), while a real division takes ``t_cycle_s``. So
    ``S = t_cycle_s / (dt · div_every_steps)``.

    Args:
        dt_s: BAOAB timestep [s].
        div_every_steps: division-check cadence in steps.
        t_cycle_s: real cell-cycle / division time [s] (default ~24 h).

    Returns:
        The dimensionless acceleration factor S (> 0).
    """
    t_sim_per_div = float(dt_s) * float(div_every_steps)
    if t_sim_per_div <= 0:
        raise ValueError("dt_s · div_every_steps must be > 0")
    return float(t_cycle_s) / t_sim_per_div


@dataclass(frozen=True, slots=True)
class DivisionRateAnchor:
    """A PHYSICALLY-ANCHORED division probability for a proliferation run (pure data).

    Derived (no free parameter, no fitting-to-PI-band) from the REAL cell-cycle
    time relative to the run's real-time-equivalent duration via the spreading-front
    acceleration factor S (``common/sim_realtime``). See :func:`division_probability`.

    Attributes:
        p_div: per-eligible-rim-cell division probability PER division-check
            (the value to feed ``GpuProliferationUpdater(p_div=…)``).
        t_real_s: run real-time equivalent = S · dt · steps [s] (the spreading
            episode's biological-equivalent length).
        t_cycle_s: real cell-cycle / division time used [s].
        div_per_cell: expected divisions a single rim cell undergoes over the run
            = t_real / T_cycle (the per-cell physical expectation p_div is set to
            reproduce; ``= n_checks · p_div``).
        expected_divisions: expected TOTAL division events over the spheroid
            = div_per_cell · n_rim_est (the spheroid-level prediction, NOT a tuned
            target). Typically ≪ 1 for a single short spreading episode.
        n_checks: number of division-checks over the run = steps // div_every.
        accel_factor: the spreading-front acceleration factor S used.
        detail: one-line human derivation note (for the run header / metrics).
    """

    p_div: float
    t_real_s: float
    t_cycle_s: float
    div_per_cell: float
    expected_divisions: float
    n_checks: int
    accel_factor: float
    detail: str

    def to_dict(self) -> dict:
        """JSON-serialisable provenance block (drop into a metrics dict)."""
        return {
            "auto_p_div": float(self.p_div),
            "div_t_real_s": float(self.t_real_s),
            "div_t_real_human": format_time(self.t_real_s),
            "cell_cycle_s": float(self.t_cycle_s),
            "cell_cycle_human": format_time(self.t_cycle_s),
            "div_per_cell": float(self.div_per_cell),
            "expected_divisions": float(self.expected_divisions),
            "div_n_checks": int(self.n_checks),
            "div_accel_factor": float(self.accel_factor),
            "div_anchor_detail": self.detail,
        }


def division_probability(steps: int, dt_s: float, *, t_cycle_h: float,
                         div_every: int, n_rim_est: int,
                         accel: float = DEFAULT_SPREADING_ACCEL
                         ) -> DivisionRateAnchor:
    """Derive the per-rim-cell division probability from PHYSICS (no free knob).

    The mechanistic anchor (PI 2026-06-11): each rim cell carries its OWN cell
    cycle of real length ``T_cycle`` (12–24 h). The run's real-time-equivalent
    duration is ``t_real = S · dt · steps`` (the spreading-front mapping, eq. 1 of
    this module). A division-check fires every ``div_every`` steps, i.e. every

        Δt_real_per_check = S · dt · div_every                                  (s)

    of real-time-equivalent. Over one check interval a single rim cell completes a
    fraction ``Δt_real_per_check / T_cycle`` of its cycle, so the PHYSICAL
    per-check division probability is

        p_div = (S · dt · div_every) / T_cycle                                  (3)

    This makes the realised divisions a PREDICTION: over ``n_checks = steps //
    div_every`` checks a rim cell's expected divisions are
    ``n_checks · p_div = t_real / T_cycle = div_per_cell`` — the cell-cycle fraction
    the spreading episode occupies — and the spheroid total is ``div_per_cell ·
    n_rim``. Because one spreading episode (seconds–hours) is ≪ one 12–24 h cycle,
    ``div_per_cell ≪ 1`` and the expected total is typically ≪ 1: proliferation
    barely moves A/A₀ over a single episode. That LOW number IS the honest physical
    answer — it is NOT chosen to hit any A/A₀ band (the b/R proliferation term in
    the PI law reflects spreading over the FULL days-long assay, not one episode).

    Args:
        steps: integrator steps the run will execute.
        dt_s: BAOAB timestep [s].
        t_cycle_h: real cell-cycle / division time [hours] (MCF7 ~18–24 h).
        div_every: division-check cadence in steps.
        n_rim_est: estimate of the number of eligible (rim) cells, for the
            expected-TOTAL-divisions prediction only (does NOT enter ``p_div``).
        accel: spreading-front acceleration factor S (default = the H.7 anchor
            6e5; the SAME S the real-time mapping uses for the spreading episode).

    Returns:
        A :class:`DivisionRateAnchor` (``p_div`` + full provenance). ``p_div`` is
        clamped to [0, 1] (it is a probability); the un-clamped derivation is kept
        in ``detail`` if the clamp ever bites (it never does at physiological
        bands — ``p_div`` is ~1e-5).
    """
    if steps < 0:
        raise ValueError(f"steps must be >= 0, got {steps}")
    if dt_s <= 0:
        raise ValueError(f"dt_s must be > 0, got {dt_s}")
    if div_every <= 0:
        raise ValueError(f"div_every must be > 0, got {div_every}")
    if t_cycle_h <= 0:
        raise ValueError(f"t_cycle_h must be > 0, got {t_cycle_h}")
    if accel <= 0:
        raise ValueError(f"accel must be > 0, got {accel}")
    t_cycle_s = float(t_cycle_h) * 3600.0
    t_real_s = float(accel) * float(dt_s) * float(steps)
    # eq. (3): per-check real-time fraction of one cell cycle.
    dt_real_per_check = float(accel) * float(dt_s) * float(div_every)
    p_raw = dt_real_per_check / t_cycle_s
    p_div = min(1.0, max(0.0, p_raw))
    n_checks = int(steps) // int(div_every)
    div_per_cell = t_real_s / t_cycle_s
    expected_divisions = div_per_cell * float(max(0, n_rim_est))
    detail = (
        f"p_div={p_raw:.2e} = S·dt·div_every/T_cycle "
        f"(S={_sci(accel)}, dt={dt_s:g}s, div_every={div_every}, "
        f"T_cycle={format_time(t_cycle_s)}); "
        f"t_real={format_time(t_real_s)} ⇒ div/cell={div_per_cell:.2e}, "
        f"expected total≈{expected_divisions:.2g} over n_rim≈{n_rim_est} "
        f"({n_checks} checks). PHYSICAL prediction, not tuned to A/A₀."
    )
    return DivisionRateAnchor(
        p_div=p_div, t_real_s=t_real_s, t_cycle_s=t_cycle_s,
        div_per_cell=div_per_cell, expected_divisions=expected_divisions,
        n_checks=n_checks, accel_factor=float(accel), detail=detail)


def map_realtime(steps: int, dt_s: float, *,
                 accel: float = DEFAULT_SPREADING_ACCEL,
                 basis: str = "spreading_front",
                 detail: str | None = None) -> RealTimeMapping:
    """Map ``(steps, dt)`` → a :class:`RealTimeMapping` (the core entry point).

    ``t_sim = dt · steps`` (raw integrator time); ``t_real = accel · t_sim``
    (real-time equivalent via the documented acceleration factor S = ``accel``).
    Pure-python, no sim cost.

    Args:
        steps: integrator steps run.
        dt_s: BAOAB timestep [s].
        accel: acceleration factor S (default = the spreading anchor, 6e5). Build
            it from :func:`spreading_accel` or :func:`proliferation_accel`.
        basis: the rate-limiting process the factor is based on
            ("spreading_front" | "cell_cycle"), recorded for provenance.
        detail: optional one-line derivation note for captions/metrics; a sane
            default is generated if omitted.

    Returns:
        A :class:`RealTimeMapping`.
    """
    if steps < 0:
        raise ValueError(f"steps must be >= 0, got {steps}")
    if dt_s <= 0:
        raise ValueError(f"dt_s must be > 0, got {dt_s}")
    if accel <= 0:
        raise ValueError(f"accel must be > 0, got {accel}")
    t_sim = float(dt_s) * float(steps)
    t_real = float(accel) * t_sim
    if detail is None:
        if basis == "cell_cycle":
            detail = (f"real-time equiv. via cell-cycle accel S={_sci(accel)} "
                      f"(one division ~ {format_time(T_CYCLE_S_DEFAULT)} real)")
        else:
            detail = (f"real-time equiv. via spreading-front accel S={_sci(accel)} "
                      f"(v_sim/v_real; v_real~{V_REAL_PROTRUSION_UM_MIN_LO}-"
                      f"{V_REAL_SPREAD_FRONT_UM_MIN_HI} µm/min). MAPPING not a "
                      f"native real-time integration.")
    return RealTimeMapping(steps=int(steps), dt_s=float(dt_s), t_sim_s=t_sim,
                           accel_factor=float(accel), t_real_s=t_real,
                           basis=str(basis), detail=str(detail))


# ---------------------------------------------------------------------------
# Human-readable time formatting
# ---------------------------------------------------------------------------
def format_time(seconds: float) -> str:
    """Format a duration in seconds to human units (µs / ms / s / min / h / d).

    Picks the largest unit in which the value is ≥ 1 (or µs for sub-ms), so a
    physical-time readout is always legible. Negative inputs keep their sign;
    non-finite inputs are passed through as their string.

    Examples:
        >>> format_time(12e-6)
        '12 µs'
        >>> format_time(0.0072)
        '7.2 ms'
        >>> format_time(7.2)
        '7.2 s'
        >>> format_time(159.6)
        '2.66 min'
        >>> format_time(7200.0)
        '2 h'
        >>> format_time(172800.0)
        '2 d'
    """
    import math

    if seconds != seconds or math.isinf(seconds):  # NaN / inf
        return str(seconds)
    sign = "-" if seconds < 0 else ""
    s = abs(float(seconds))
    if s == 0.0:
        return "0 s"
    minute, hour, day = 60.0, 3600.0, 86400.0
    if s < 1.0e-3:
        return f"{sign}{_fmt_num(s * 1.0e6)} µs"
    if s < 1.0:
        return f"{sign}{_fmt_num(s * 1.0e3)} ms"
    if s < minute:
        return f"{sign}{_fmt_num(s)} s"
    if s < hour:
        return f"{sign}{_fmt_num(s / minute)} min"
    if s < day:
        return f"{sign}{_fmt_num(s / hour)} h"
    return f"{sign}{_fmt_num(s / day)} d"


def _fmt_num(x: float) -> str:
    """Compact number: integer if whole, else up to 3 significant-ish digits."""
    if x == int(x):
        return str(int(x))
    if x >= 100:
        return f"{x:.0f}"
    if x >= 10:
        return f"{x:.1f}"
    return f"{x:.2f}".rstrip("0").rstrip(".")


def _sci(x: float) -> str:
    """Render a factor like 6e5 as ``6×10⁵`` (superscript) for titles."""
    import math

    if x <= 0 or x != x or math.isinf(x):
        return str(x)
    exp = int(math.floor(math.log10(x)))
    mant = x / (10.0 ** exp)
    sup = str(exp).translate(str.maketrans("0123456789-", "⁰¹²"
                                           "³⁴⁵⁶⁷"
                                           "⁸⁹⁻"))
    if abs(mant - 1.0) < 1e-9:
        return f"10{sup}"
    return f"{_fmt_num(mant)}×10{sup}"
