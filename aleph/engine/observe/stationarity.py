"""Has a driven steady state actually been reached, and what is the error bar on its mean?

WHY THIS MODULE EXISTS.  Nothing in this repository could answer either question.  Every gate scored a
static residual or a fixed step count, which is the shape you use when you believe the system relaxes to
rest.  It does not: a cortex's tension is maintained by myosin turnover, so the object to measure is a
**non-equilibrium steady state** — statistically stationary, never still.

The cost of not having this is on the record.  GATE B reported "cortical tension emerges 0 → 3.72 pN/µm"
from 30 steps at ``dt_phys`` 0.01 — **0.3 s of physical time against a 2.5 s bound-head lifetime**
(``k_off0`` = 0.4/s, Nagy 2013).  It observed 12% of one turnover and called it a measurement, and no
check existed that could have said otherwise.  :class:`StationarityVerdict.TOO_SHORT` exists so that
exact run returns a refusal rather than a number.

THREE THINGS IT DOES THAT A NAIVE MEAN DOES NOT.

1. **It measures the correlation time instead of assuming it.**  In a driven state the collective
   relaxation can be far slower than the single-head lifetime that drives it, so ``1/k_off`` is used
   ONLY to set the shortest series length worth judging — never as the correlation time itself.
2. **It reports the error on the mean using the number of INDEPENDENT samples.**  Consecutive steps of a
   correlated series are not independent; ``sem = σ·sqrt(2·τ_int/T)``, not ``σ/sqrt(N)``.  With
   ``dt_phys`` = 0.01 s and a correlation time of order a second, those differ by more than an order of
   magnitude, and the naive form is the standard way to publish an error bar that is far too small.
3. **It chooses where the transient ends by a rule, not by eye** — the discard point that maximises the
   number of independent samples in what remains (Chodera 2016).  This project has retracted enough
   numbers that were read off a plot; the discard point is now an output, with its own value recorded.

WHAT IT DELIBERATELY DOES NOT DO.  It does not assume the fluctuation–dissipation theorem.  In a driven
steady state FDT does not hold, so the fluctuation size cannot be predicted from temperature and is
measured.  It also takes no position on whether a stationary value is *correct* — that is the gate's
band, a separate question from whether a mean is meaningful at all.

engine units: time in SECONDS of physical time, never solver iterations.

Sanity Gate:
    * dimensional: ``dt`` and every returned time are [s]; ``tau_int`` is [s]; ``n_eff`` is dimensionless.
      A caller passing iteration index as time gets a τ in "iterations", which is why ``dt`` is required
      and has no default.
    * boundary: a series shorter than ``min_windows`` correlation times returns ``TOO_SHORT`` with no
      mean — the GATE B case.  A constant series has zero variance and returns τ_int = dt with a note,
      rather than dividing by zero.
    * conservation/invariant: ``n_eff <= n_samples`` always, with equality only for an uncorrelated
      series; the automatic window is the Sokal criterion, so τ_int is not summed into its own noise.
    * numerical: the autocorrelation is computed by FFT, so the cost is ``N log N`` and a 2,500-step
      series is free; the mean is subtracted before correlating, since not doing so measures the mean
      rather than the fluctuation.
    * sign-sense: drift is signed and reported signed — a tension falling toward its steady state and one
      rising toward it are different physics, and ``abs`` would hide which happened.
    * measurement-protocol: ``min_windows`` and ``drift_sigma`` are the caller's PRE-DECLARED contract,
      recorded in the report next to the verdict, because choosing them after seeing the series is how a
      transient gets certified as a steady state.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any

import numpy as np

__all__ = [
    "StationarityVerdict",
    "StationarityReport",
    "integrated_autocorrelation_time",
    "equilibration_point",
    "assess_stationarity",
    "assess_observables",
]


class StationarityVerdict(StrEnum):
    """Whether a time series supports quoting a steady-state mean.

    * ``STATIONARY`` — the trend over the analysed window is smaller than the series' own fluctuation,
      and the window is long enough relative to the measured correlation time for that to mean something.
    * ``DRIFTING`` — the series is long enough to judge and is still going somewhere.  A mean here is the
      average of a transient, which is what the withdrawn tension numbers were.
    * ``TOO_SHORT`` — not enough physical time to distinguish the two.  Distinct from ``DRIFTING`` on
      purpose: "we looked and it moved" and "we did not look long enough" are different findings, and
      collapsing them is what let 0.3 s of a 2.5 s process be reported as a measurement.
    """

    STATIONARY = "STATIONARY"
    DRIFTING = "DRIFTING"
    TOO_SHORT = "TOO_SHORT"


@dataclass(frozen=True, slots=True)
class StationarityReport:
    """One observable's verdict, with everything needed to re-judge it without the raw series."""

    verdict: StationarityVerdict
    observable: str
    #: Steady-state mean over the analysed window — ``None`` unless the verdict is ``STATIONARY``.
    mean: float | None
    #: Standard error of that mean, computed from INDEPENDENT samples, not from ``n_samples``.
    sem: float | None
    std: float
    tau_int_s: float
    n_eff: float
    #: Physical time discarded as transient, and the rule that chose it.
    equilibration_time_s: float
    analysed_time_s: float
    #: Signed least-squares slope over the analysed window [observable units per second].
    drift_per_s: float
    #: ``|drift over the window| / std`` — what the verdict is thresholded on.
    drift_over_std: float
    contract: dict[str, Any] = field(default_factory=dict)
    notes: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        """JSON-able form for an :func:`observation_artifact` measurements block."""
        return {
            "verdict": self.verdict.value, "observable": self.observable,
            "mean": self.mean, "sem": self.sem, "std": self.std,
            "tau_int_s": self.tau_int_s, "n_eff": self.n_eff,
            "equilibration_time_s": self.equilibration_time_s,
            "analysed_time_s": self.analysed_time_s,
            "drift_per_s": self.drift_per_s, "drift_over_std": self.drift_over_std,
            "contract": dict(self.contract), "notes": dict(self.notes),
        }


def integrated_autocorrelation_time(
    series: np.ndarray, dt: float, *, c: float = 5.0
) -> tuple[float, float, str]:
    """Return ``(tau_int_s, n_eff, note)`` for one series.

    ``tau_int = dt · (1/2 + Σ_{k=1..M} ρ_k)`` with ρ the normalised autocorrelation and ``M`` chosen by
    the Sokal automatic window: the smallest ``M`` with ``M >= c · tau_int(M)``.  The window is the whole
    point — the tail of ρ is noise of order ``1/sqrt(N)`` per lag, so summing every lag adds a random
    walk to the answer and typically inflates τ by more than the signal.

    Args:
        series: Observable samples, evenly spaced in time.
        dt: Physical time between samples [s].  Required and unitful: passing an iteration index here
            silently returns a correlation time in iterations, which is exactly the confusion this whole
            module exists to prevent.
        c: Sokal window constant.  6 is the usual choice for MD; 5 is used here because the series this
            judges are short by MD standards and a larger ``c`` refuses to close the window at all.

    Returns:
        ``tau_int`` in seconds (floored at ``dt/2``, the value for a perfectly uncorrelated series),
        the effective independent-sample count ``n_eff = T / (2·tau_int)``, and a note naming what
        happened when the estimate was degenerate.

    Raises:
        ValueError: If ``dt`` is not positive-finite or the series has fewer than 4 samples.
    """
    if not (dt > 0.0 and math.isfinite(dt)):
        raise ValueError(f"dt must be positive-finite seconds; got {dt!r}")
    x = np.asarray(series, dtype=np.float64).ravel()
    if x.size < 4:
        raise ValueError(f"need >= 4 samples to estimate a correlation time; got {x.size}")

    n = x.size
    total_time = n * dt
    var = float(np.var(x))
    if var <= 0.0 or not math.isfinite(var):
        return dt * 0.5, float(n), "series has zero variance; tau_int floored at dt/2"

    # FFT autocorrelation of the MEAN-SUBTRACTED series. Without the subtraction this measures the mean.
    y = x - x.mean()
    size = 1 << (2 * n - 1).bit_length()
    spectrum = np.fft.rfft(y, size)
    acf = np.fft.irfft(spectrum * np.conjugate(spectrum), size)[:n].real
    acf /= acf[0]

    tau = 0.5
    window = n - 1
    for m in range(1, n):
        tau = 0.5 + float(np.sum(acf[1 : m + 1]))
        if m >= c * tau:
            window = m
            break
    else:
        note = (f"Sokal window never closed within {n} samples: the series is short relative to its own "
                "correlation time, so tau_int is a LOWER bound and n_eff an upper one")
        tau = max(tau, 0.5)
        return tau * dt, max(1.0, total_time / (2.0 * tau * dt)), note

    tau = max(tau, 0.5)
    tau_s = tau * dt
    return tau_s, max(1.0, total_time / (2.0 * tau_s)), f"Sokal window closed at lag {window} (c={c})"


def equilibration_point(series: np.ndarray, dt: float, *, n_candidates: int = 40) -> tuple[int, float]:
    """Return ``(start_index, n_eff)`` for the discard point that maximises independent samples.

    Chodera's rule: sweep candidate truncations, and keep the one whose REMAINING series has the most
    effective samples.  Discarding too little leaves the transient in and inflates the correlation time;
    discarding too much throws away signal.  The maximum trades those off without anyone looking at a
    plot, which matters here specifically — several retracted numbers in this project were read off one.

    Args:
        series: Observable samples.
        dt: Physical time between samples [s].
        n_candidates: How many truncation points to try, spread over the first half of the series.

    Returns:
        The index to start from, and the ``n_eff`` achieved there.
    """
    x = np.asarray(series, dtype=np.float64).ravel()
    best_index, best_n_eff = 0, -1.0
    limit = max(1, x.size // 2)
    step = max(1, limit // max(1, n_candidates))
    for start in range(0, limit, step):
        tail = x[start:]
        if tail.size < 4:
            break
        try:
            _, n_eff, _ = integrated_autocorrelation_time(tail, dt)
        except ValueError:
            continue
        if n_eff > best_n_eff:
            best_index, best_n_eff = start, n_eff
    return best_index, max(best_n_eff, 0.0)


def window_opens_on_a_transient(series: np.ndarray, start_index: int) -> tuple[bool, float]:
    """Return ``(is_contaminated, ratio)`` for the window :func:`equilibration_point` chose.

    WHY THIS EXISTS.  :func:`equilibration_point` maximises ``n_eff = T / (2·tau_int)`` on the assumption,
    stated in its own docstring, that "discarding too little leaves the transient in and INFLATES the
    correlation time".  On 2026-07-29 that failed on a real run: a series rising 0.8 -> 3.4 in 0.5 s and
    then flat for 7 s got ``start_index = 0``.  ``T`` is maximal there and the ACF estimator did not
    penalise the trend enough to compensate, so the "steady-state mean" was taken ACROSS THE RISE — the
    exact failure this module exists to prevent.

    It also self-reinforced.  The drift verdict is ``|drift| < drift_sigma · std``, and keeping the
    transient inflated ``std`` twentyfold (0.0106 -> 0.2538) — so a LARGER transient makes the drift test
    EASIER to pass.  That run was reported STATIONARY at 0.94 sigma.

    THE STATISTIC, and why it is this one.  A transient is a property of a SEGMENT, not of a point: the
    first draft tested the first SAMPLE against the window's scatter and rejected all three real runs,
    because a cut placed at the knee necessarily leaves its first sample in the tail.  This compares the
    MEAN of the leading tenth against the mean of the trailing tenth, and normalises by the scatter of the
    TRAILING tenth — which is transient-free by construction, so the inflation loop above cannot weaken the
    test the way it weakens the drift test.  The threshold is one standard deviation: a definition, not a
    number chosen after seeing the data.  On the three runs it reads 0.0 / 0.8 / 33.5.

    Args:
        series: The full observable series.
        start_index: The discard point :func:`equilibration_point` selected.

    Returns:
        ``(True, ratio)`` when the retained window still opens on a rise, else ``(False, ratio)``.
    """
    x = np.asarray(series, dtype=np.float64).ravel()
    tail = x[start_index:]
    if tail.size < 20:                       # too short to split into tenths and say anything
        return (False, 0.0)
    k = max(2, tail.size // 10)
    lead = float(np.mean(tail[:k]))
    trail_seg = tail[-k:]
    trail = float(np.mean(trail_seg))
    spread = float(np.std(trail_seg, ddof=1))
    if not np.isfinite(spread) or spread <= 0.0:
        return (False, 0.0)
    ratio = abs(lead - trail) / spread
    return (ratio > 1.0, ratio)


def assess_stationarity(
    series: np.ndarray,
    dt: float,
    *,
    observable: str,
    driving_correlation_time_s: float,
    min_windows: float = 5.0,
    min_tau_windows: float = 50.0,
    drift_sigma: float = 1.0,
) -> StationarityReport:
    """Judge one observable, returning a mean ONLY if a steady state was actually reached.

    The two thresholds are the caller's contract and must be declared before the run; they are copied
    into the report so a later reader judges the verdict against the criterion that produced it.

    Args:
        series: Observable samples, evenly spaced.
        dt: Physical time between samples [s].
        observable: Name, recorded in the report.
        driving_correlation_time_s: The physical time the DRIVING process sets — here the NMII bound-head
            lifetime ``1/k_off0`` ≈ 2.5 s.  Used only to reject a series too short to judge; the
            observable's own correlation time is measured, because collective relaxation in a driven
            state is routinely much slower than the process driving it.
        min_windows: Shortest analysable span, in units of ``driving_correlation_time_s``, that may
            return anything other than ``TOO_SHORT``.  This is an EXTERNAL bound: it knows nothing about
            the observable's own relaxation.
        min_tau_windows: The SELF-REFERENTIAL bound, in units of the MEASURED ``tau_int``.  Equivalent to
            ``n_eff >= min_tau_windows / 2`` — at least 25 independent samples at the default 50.

            ⚠ **Adopted 2026-08-20 by PI decision** ("일단 B로 진행"), from
            :mod:`aleph.observe.stationarity`, whose docstring names the external-bound-only design this
            function previously used as a failure mode: *"the driving time is a property of the forcing,
            and the collective relaxation of the observable can be orders of magnitude slower. A series
            can clear that bar and still be hopeless for its own correlation time."*  Measured before the
            change: at the NMII driving time of 2.5 s the external bound requires 12.5 s, while 50 x the
            measured tau ranges from 0.08 s to 29 s over AR(1) rho 0.5-0.999 — so neither bound dominates
            and the two are not the same knob.  Both are kept, exactly as that module recommends ("a
            caller who genuinely has an external bound can impose it as well"), and both are recorded in
            the returned contract.
        drift_sigma: Stationary requires ``|drift across the window| < drift_sigma · std``.  **LEFT AT
            1.0, and NOT raised to the adopted module's 3.0.**  The PI's decision was taken on a
            comparison this session got wrong — ``min_windows`` 5 and ``min_tau_windows`` 50 were
            presented as the same knob, and they are not.  On this knob the adopted contract is LOOSER,
            and moving it measurably flips cases: at rho=0.9 with a linear drift the window scores
            ``|drift|/std = 1.00``, which 3.0 admits and 1.0 refuses.  Raising it moved the dominant
            disagreement with :mod:`aleph.observe.stationarity` from 74 cases to 89 — i.e. loosening
            here made the two modules agree LESS, not more.  A gate threshold moved on a mistaken premise
            is what the charter forbids, so this half is held and surfaced rather than applied.

    Returns:
        The :class:`StationarityReport`.

    Raises:
        ValueError: If ``dt`` or ``driving_correlation_time_s`` is not positive-finite.
    """
    if not (dt > 0.0 and math.isfinite(dt)):
        raise ValueError(f"dt must be positive-finite seconds; got {dt!r}")
    if not (driving_correlation_time_s > 0.0 and math.isfinite(driving_correlation_time_s)):
        raise ValueError(
            f"driving_correlation_time_s must be positive-finite; got {driving_correlation_time_s!r}"
        )

    x = np.asarray(series, dtype=np.float64).ravel()
    contract = {
        "min_windows": float(min_windows), "min_tau_windows": float(min_tau_windows),
        "drift_sigma": float(drift_sigma),
        "driving_correlation_time_s": float(driving_correlation_time_s),
        "declared": "before the run; recorded here so the verdict is judged against its own criterion",
    }
    required_s = float(min_windows) * float(driving_correlation_time_s)

    if x.size < 4 or x.size * dt < required_s:
        return StationarityReport(
            verdict=StationarityVerdict.TOO_SHORT, observable=observable, mean=None, sem=None,
            std=float(np.std(x)) if x.size else 0.0, tau_int_s=float("nan"), n_eff=0.0,
            equilibration_time_s=0.0, analysed_time_s=float(x.size * dt),
            drift_per_s=float("nan"), drift_over_std=float("nan"), contract=contract,
            notes={"reason": (
                f"{x.size * dt:.3g} s of physical time against a required {required_s:.3g} s "
                f"({min_windows:g} x the {driving_correlation_time_s:.3g} s driving correlation time). "
                "This is the GATE B case: a transient cannot be told from a steady state here, so no "
                "mean is returned")},
        )

    start, _ = equilibration_point(x, dt)
    contaminated, open_ratio = window_opens_on_a_transient(x, start)
    if contaminated:
        # The detector kept a window that still opens on a transient — see that function for the run this
        # was written against. Report it rather than silently choosing a different cut: a truncation the
        # estimator did not pick is one a human picked after seeing the answer, which is the thing this
        # module is for preventing. No mean is returned, so nothing downstream can quote the contaminated
        # average that would otherwise have been produced.
        kept = x[start:]
        # Drift is still reported SIGNED over the retained window. The module's contract is that a series
        # rising toward its steady state and one falling toward it are different physics, so `abs` (or a
        # nan) would hide which happened — and a caller who only learns "DRIFTING" cannot tell whether to
        # run longer or to look for an instability.
        kept_std = float(np.std(kept))
        kept_times = np.arange(kept.size, dtype=np.float64) * dt
        kept_slope = float(np.polyfit(kept_times, kept, 1)[0]) if kept.size >= 2 else float("nan")
        kept_ratio = (abs(kept_slope * kept.size * dt) / kept_std) if kept_std > 0.0 else 0.0
        return StationarityReport(
            verdict=StationarityVerdict.DRIFTING, observable=observable, mean=None, sem=None,
            std=kept_std, tau_int_s=float("nan"), n_eff=0.0,
            equilibration_time_s=float(start * dt), analysed_time_s=float(kept.size * dt),
            drift_per_s=kept_slope, drift_over_std=kept_ratio, contract=contract,
            notes={"reason": (
                f"the equilibration search kept t >= {start * dt:.3g} s, but that window still OPENS on a "
                f"transient: its leading tenth sits {open_ratio:.1f} standard deviations below its "
                "trailing tenth. Averaging across a rise is what this detector exists to prevent, and "
                "keeping one also INFLATES std, which makes the drift test EASIER to pass — the verdict "
                "would have been self-reinforcing. Run longer, or sample the transient more finely so the "
                "search has candidates inside it")},
        )
    tail = x[start:]
    analysed_s = float(tail.size * dt)
    if analysed_s < required_s:
        # The discard ate the window. Before conceding TOO_SHORT, ask the FULL series: if it is long
        # enough to judge and is plainly still going somewhere, "it is still drifting" is a conclusion we
        # are entitled to draw, and reporting TOO_SHORT instead would hide a real finding behind a
        # bookkeeping rule. Only a series that is neither clearly drifting nor long enough to call
        # stationary is genuinely inconclusive.
        full_times = np.arange(x.size, dtype=np.float64) * dt
        full_std = float(np.std(x))
        full_slope = float(np.polyfit(full_times, x, 1)[0])
        full_ratio = abs(full_slope * x.size * dt) / full_std if full_std > 0.0 else 0.0
        if full_ratio >= float(drift_sigma):
            return StationarityReport(
                verdict=StationarityVerdict.DRIFTING, observable=observable, mean=None, sem=None,
                std=full_std, tau_int_s=float("nan"), n_eff=0.0,
                equilibration_time_s=float(start * dt), analysed_time_s=float(x.size * dt),
                drift_per_s=full_slope, drift_over_std=full_ratio, contract=contract,
                notes={"reason": (
                    "judged on the FULL series: no truncation leaves a long enough window, and the "
                    f"whole span drifts by {full_ratio:.2f}x its own fluctuation, so this is still a "
                    "transient rather than an inconclusive one")},
            )
        return StationarityReport(
            verdict=StationarityVerdict.TOO_SHORT, observable=observable, mean=None, sem=None,
            std=float(np.std(tail)), tau_int_s=float("nan"), n_eff=0.0,
            equilibration_time_s=float(start * dt), analysed_time_s=analysed_s,
            drift_per_s=float("nan"), drift_over_std=float("nan"), contract=contract,
            notes={"reason": (
                f"after discarding {start * dt:.3g} s of transient only {analysed_s:.3g} s remain, "
                f"under the required {required_s:.3g} s. Run longer: the budget must cover the "
                "transient PLUS the window, not the window alone")},
        )

    tau_s, n_eff, tau_note = integrated_autocorrelation_time(tail, dt)

    # The self-referential refusal. Everything above tested the window against a time scale supplied from
    # OUTSIDE — the driving process — which a series can clear while still holding almost no independent
    # samples of the thing being measured. This asks the series about itself.
    if analysed_s < float(min_tau_windows) * tau_s:
        return StationarityReport(
            verdict=StationarityVerdict.TOO_SHORT, observable=observable, mean=None, sem=None,
            std=float(np.std(tail)), tau_int_s=tau_s, n_eff=n_eff,
            equilibration_time_s=float(start * dt), analysed_time_s=analysed_s,
            drift_per_s=float("nan"), drift_over_std=float("nan"), contract=contract,
            notes={"reason": (
                f"{analysed_s:.3g} s of analysed window against a required "
                f"{float(min_tau_windows) * tau_s:.3g} s ({min_tau_windows:g} x the MEASURED "
                f"tau_int of {tau_s:.3g} s), i.e. n_eff = {n_eff:.1f} against a required "
                f"{float(min_tau_windows) / 2.0:.0f} independent samples. The external driving-time "
                "bound was cleared and this one was not, which is the case that bound cannot see"),
                "tau_estimate": tau_note,
                "bias_caveat": (
                    "tau_int measured on a short series is BIASED LOW, so this refusal is if anything "
                    "too permissive — see aleph/observe/stationarity.py for the measured bias table"),
            },
        )

    std = float(np.std(tail))
    times = np.arange(tail.size, dtype=np.float64) * dt
    slope = float(np.polyfit(times, tail, 1)[0])
    drift_across_window = slope * analysed_s
    drift_over_std = abs(drift_across_window) / std if std > 0.0 else 0.0

    stationary = drift_over_std < float(drift_sigma)
    mean = float(np.mean(tail)) if stationary else None
    sem = (std * math.sqrt(2.0 * tau_s / analysed_s)) if stationary else None

    return StationarityReport(
        verdict=StationarityVerdict.STATIONARY if stationary else StationarityVerdict.DRIFTING,
        observable=observable, mean=mean, sem=sem, std=std, tau_int_s=tau_s, n_eff=n_eff,
        equilibration_time_s=float(start * dt), analysed_time_s=analysed_s,
        drift_per_s=slope, drift_over_std=drift_over_std, contract=contract,
        notes={
            "tau_estimate": tau_note,
            "sem_basis": (
                "sigma*sqrt(2*tau_int/T) — the error on a correlated mean. The naive sigma/sqrt(N) would "
                f"be {std / math.sqrt(tail.size):.4g} here, i.e. "
                f"{(std * math.sqrt(2.0 * tau_s / analysed_s)) / (std / math.sqrt(tail.size)):.3g}x too "
                "small" if std > 0.0 else "zero variance"),
            "fdt": ("not assumed — in a driven steady state the fluctuation size is not set by "
                    "temperature, so std is measured rather than predicted"),
        },
    )


def assess_observables(
    trajectory: list[dict],
    sample_dt_s: float,
    *,
    observables: tuple[str, ...],
    driving_correlation_time_s: float,
    min_windows: float = 5.0,
    drift_sigma: float = 1.0,
) -> dict[str, dict[str, Any]]:
    """Judge every named column of a recorded trajectory and return JSON-able reports.

    WHY THIS IS SHARED CODE.  A driver that writes a trajectory has to do exactly this, and the first
    driver to need it grew its own copy — which is the pattern that let a wrong bound-head lifetime
    survive in one driver while the others could not even express the question.  A second driver should
    add this in three lines, not ninety.

    ``sample_dt_s`` is the SAMPLE spacing, not ``dt_phys``.  With a telemetry stride the two differ, and
    handing the physical step to an estimator that is fed every stride-th sample understates the
    correlation time by exactly that stride — which would make ``n_eff`` and the error bar too good in
    the dangerous direction.

    Args:
        trajectory: Recorded rows; each must carry every name in ``observables``.
        sample_dt_s: Physical time between consecutive rows [s].
        observables: Column names to judge, in report order.
        driving_correlation_time_s: The driving time the contract is expressed in [s] — for a motor
            lane, :func:`~aleph.components.motor.bell_kinetics_analytic.bound_state_lifetime`'s conservative
            value, never a literal.
        min_windows: Contract — shortest analysable span in units of ``driving_correlation_time_s``.
        drift_sigma: Contract — stationary requires drift across the window below this many sigma.

    Returns:
        ``{observable: report}``.  A column with too few samples to judge is reported ``TOO_SHORT``
        rather than skipped, because a missing verdict reads as an absent problem.

    Raises:
        KeyError: If a row lacks one of ``observables`` — a silently dropped column would be a silently
            unjudged one.
    """
    reports: dict[str, dict[str, Any]] = {}
    for name in observables:
        series = np.array([row[name] for row in trajectory], dtype=np.float64)
        if series.size < 4:
            reports[name] = {
                "verdict": StationarityVerdict.TOO_SHORT.value, "observable": name,
                "mean": None, "sem": None, "n_samples": int(series.size),
                "notes": {"reason": f"{series.size} samples; the estimator needs >= 4"},
            }
            continue
        report = assess_stationarity(
            series, sample_dt_s, observable=name,
            driving_correlation_time_s=driving_correlation_time_s,
            min_windows=min_windows, drift_sigma=drift_sigma)
        reports[name] = {**report.as_dict(), "n_samples": int(series.size)}
    return reports
