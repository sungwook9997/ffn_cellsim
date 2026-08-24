"""The stationarity detector must be able to say TOO_SHORT, DRIFTING and STATIONARY — and be right.

WHY THESE TESTS EXIST.  ``observe/stationarity.py`` landed on 2026-07-28 as the guard against reporting
the average of a transient, and it landed with **no tests at all**.  A detector nobody has shown can
distinguish its three verdicts is not a guard; it is a second opinion of unknown quality standing in
front of a gate.  The cases below are chosen so that each verdict is forced by a series whose right
answer is known analytically rather than by inspection:

* white noise has ``tau_int = dt/2`` exactly, so a detector that inflates tau fails here;
* AR(1) has ``tau_int = dt(1+phi)/(2(1-phi))`` in closed form, so the Sokal window is checked against a
  number rather than against plausibility;
* a ramp is unambiguously DRIFTING at any length, which is the verdict the withdrawn tension numbers
  should have received;
* the GATE B case — 0.3 s of physical time against a 2.5 s bound-head lifetime — must come back
  TOO_SHORT and must NOT come back DRIFTING, because "we did not look long enough" and "we looked and
  it moved" are different findings and collapsing them is the specific error being guarded against.

Runs on CPU: this module is pure numpy analysis of a recorded series and touches no device.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from aleph.engine.observe.stationarity import (
    StationarityReport,
    StationarityVerdict,
    assess_stationarity,
    equilibration_point,
    integrated_autocorrelation_time,
)

#: The NMII bound-head lifetime the cortex runs are judged against [s].  Kept explicit here rather than
#: imported so a change to the engine constant shows up as a test edit, not as a silently moved goalpost.
DRIVING_TAU_S = 1.0 / 0.4


def _ar1(n: int, phi: float, *, seed: int, sigma: float = 1.0) -> np.ndarray:
    """Return an AR(1) series ``x_{k+1} = phi*x_k + eps``, started from its stationary distribution.

    Starting from the stationary variance ``sigma^2/(1-phi^2)`` matters: started at zero the series has
    a burn-in that the correlation-time estimator would legitimately read as a transient, which is a
    different test from the one intended here.
    """
    rng = np.random.default_rng(seed)
    x = np.empty(n, dtype=np.float64)
    x[0] = rng.normal(0.0, sigma / math.sqrt(1.0 - phi * phi))
    for k in range(1, n):
        x[k] = phi * x[k - 1] + rng.normal(0.0, sigma)
    return x


class TestIntegratedAutocorrelationTime:
    """tau_int must be measured, and it must be right on series whose tau is known in closed form."""

    def test_white_noise_gives_tau_of_half_a_step(self) -> None:
        """An uncorrelated series has ``tau_int = dt/2`` and ``n_eff = N``; anything larger is inflation."""
        dt = 0.01
        x = np.random.default_rng(0).normal(size=20_000)
        tau_s, n_eff, _ = integrated_autocorrelation_time(x, dt)
        assert tau_s == pytest.approx(dt * 0.5, rel=0.25)
        assert n_eff == pytest.approx(x.size, rel=0.25)

    @pytest.mark.parametrize("phi", [0.8, 0.9, 0.95])
    def test_ar1_recovers_the_closed_form_tau(self, phi: float) -> None:
        """``tau_int = dt*(1+phi)/(2*(1-phi))`` for AR(1) — the estimator is checked against algebra."""
        dt = 0.01
        expected = dt * (1.0 + phi) / (2.0 * (1.0 - phi))
        x = _ar1(200_000, phi, seed=1)
        tau_s, n_eff, _ = integrated_autocorrelation_time(x, dt)
        assert tau_s == pytest.approx(expected, rel=0.15)
        # n_eff = T/(2*tau) and is strictly fewer than the sample count for a correlated series.
        assert n_eff == pytest.approx(x.size * dt / (2.0 * expected), rel=0.2)
        assert n_eff < x.size

    def test_constant_series_does_not_divide_by_zero(self) -> None:
        """Zero variance is a boundary case, not an exception: tau floors at dt/2 with a note."""
        tau_s, n_eff, note = integrated_autocorrelation_time(np.full(64, 3.5), 0.01)
        assert tau_s == pytest.approx(0.005)
        assert n_eff == 64.0
        assert "zero variance" in note

    def test_dt_must_be_unitful_and_positive(self) -> None:
        """``dt`` has no default on purpose — an iteration index here silently returns tau in iterations."""
        x = np.random.default_rng(0).normal(size=64)
        with pytest.raises(ValueError, match="dt must be positive-finite"):
            integrated_autocorrelation_time(x, 0.0)
        with pytest.raises(ValueError, match="need >= 4 samples"):
            integrated_autocorrelation_time(np.zeros(3), 0.01)

    def test_the_window_always_closes_so_the_lower_bound_note_is_unreachable(self) -> None:
        """The "Sokal window never closed" branch is DEAD CODE, and the reason is an identity.

        The autocorrelation is taken of the MEAN-SUBTRACTED series, and for any zero-mean ``y`` the full
        two-sided sum of the unnormalised autocorrelation is ``(sum y)^2 = 0``.  Normalising by
        ``acf[0]`` therefore forces ``sum_{k>=1} rho_k = -1/2`` EXACTLY, so ``tau(n-1) = 0`` and the
        criterion ``M >= c*tau(M)`` is satisfied at the last lag whatever the series looks like.

        This matters beyond tidiness: the module's only mechanism for SAYING a series is too short for
        its own correlation time can never fire, so a short series gets a confident-looking tau instead
        of a labelled bound — see the next test for how wrong that tau can be.
        """
        for n, phi in ((8, 0.995), (24, 0.98), (40, 0.9), (80, 0.995)):
            _, _, note = integrated_autocorrelation_time(_ar1(n, phi, seed=2), 0.01)
            assert "never closed" not in note, f"branch reached at n={n}, phi={phi} — update this test"

        y = np.random.default_rng(2).normal(size=512)
        y = y - y.mean()
        size = 1 << (2 * y.size - 1).bit_length()
        spectrum = np.fft.rfft(y, size)
        acf = np.fft.irfft(spectrum * np.conjugate(spectrum), size)[: y.size].real
        assert float(np.sum(acf[1:] / acf[0])) == pytest.approx(-0.5, abs=1e-9)

    @pytest.mark.parametrize(
        "n_samples, phi, worst_ratio", [(60, 0.98, 0.20), (200, 0.98, 0.45), (1_000, 0.98, 0.80)]
    )
    def test_tau_is_underestimated_on_a_short_series_and_the_bias_is_optimistic(
        self, n_samples: int, phi: float, worst_ratio: float
    ) -> None:
        """A short series returns a tau far BELOW the truth, so ``n_eff`` and the error bar are too good.

        Measured here: at 60 samples of an AR(1) with tau_int = 49.5 steps the estimator returns under a
        fifth of it.  The bias has a direction, and it is the dangerous one — too small a tau means too
        many "independent" samples and too tight a ``sem``, which is a series being certified as settled
        on the strength of not having been watched long enough to know otherwise.

        Consequence for the run contract: a criterion expressed only in the DRIVING correlation time
        does not protect against this, because the accuracy of tau is set by the OBSERVABLE's own
        correlation time.  A window must also be long against the tau it measures — roughly 50x by these
        numbers — for the returned tau to be trustworthy.
        """
        dt = 0.01
        true_tau = dt * (1.0 + phi) / (2.0 * (1.0 - phi))
        measured = np.mean([integrated_autocorrelation_time(_ar1(n_samples, phi, seed=s), dt)[0]
                            for s in range(8)])
        assert measured < worst_ratio * true_tau, (
            f"{n_samples} samples returned tau = {measured:.4g} s vs a true {true_tau:.4g} s; if this "
            "now passes, the estimator changed and the run contract should be revisited")


class TestEquilibrationPoint:
    """The transient discard point is chosen by a rule, because reading it off a plot is what failed."""

    def test_a_decaying_transient_is_discarded(self) -> None:
        """An exponential approach to a plateau must give a strictly positive discard index."""
        dt = 0.01
        t = np.arange(4_000) * dt
        rng = np.random.default_rng(3)
        x = 5.0 - 4.0 * np.exp(-t / 0.5) + rng.normal(0.0, 0.02, size=t.size)
        start, n_eff = equilibration_point(x, dt)
        assert start > 0, "the transient was not discarded at all"
        assert start * dt < t[-1] / 2.0, "the rule may only search the first half"
        assert n_eff > 0.0

    def test_a_stationary_series_keeps_almost_everything(self) -> None:
        """With no transient the rule should not throw signal away; discarding is not free."""
        dt = 0.01
        x = np.random.default_rng(4).normal(size=8_000)
        start, _ = equilibration_point(x, dt)
        assert start * dt < 0.25 * x.size * dt


class TestAssessStationarity:
    """The three verdicts, each forced by a series whose correct answer is not a matter of opinion."""

    def test_gate_b_case_is_too_short_and_returns_no_mean(self) -> None:
        """0.3 s against a 2.5 s driving time: the exact run this module was written to refuse.

        GATE B reported "cortical tension emerges 0 -> 3.72 pN/um" from 30 steps at dt_phys 0.01 — 12%
        of one bound-head turnover.  The verdict must be TOO_SHORT, the mean must be ``None``, and the
        reason must quote both times so the refusal is self-explaining in the artifact.
        """
        dt = 0.01
        gamma = np.linspace(0.0, 3.72, 30)          # the emergence trace, 0.3 s of physical time
        report = assess_stationarity(gamma, dt, observable="gamma_total_pn_per_um",
                                     driving_correlation_time_s=DRIVING_TAU_S)
        assert report.verdict is StationarityVerdict.TOO_SHORT
        assert report.mean is None and report.sem is None
        assert report.analysed_time_s == pytest.approx(0.30)
        assert "12.5 s" in report.notes["reason"] or "12.5" in report.notes["reason"]

    def test_a_ramp_long_enough_to_judge_is_drifting_not_too_short(self) -> None:
        """A long series that is plainly still going somewhere is a CONCLUSION, not an inconclusive."""
        dt = 0.01
        x = np.linspace(0.0, 10.0, 4_000) + np.random.default_rng(5).normal(0.0, 0.05, size=4_000)
        report = assess_stationarity(x, dt, observable="ramp",
                                     driving_correlation_time_s=DRIVING_TAU_S)
        assert report.verdict is StationarityVerdict.DRIFTING
        assert report.mean is None, "a mean here is the average of a transient"
        assert report.drift_per_s > 0.0, "drift is reported SIGNED; a rising series must read positive"

    def test_stationary_noise_returns_a_mean_with_a_correlated_error_bar(self) -> None:
        """The mean is returned, and its error bar uses INDEPENDENT samples, not the sample count."""
        dt = 0.01
        phi = 0.9
        x = 4.0 + _ar1(40_000, phi, seed=6, sigma=0.1)
        report = assess_stationarity(x, dt, observable="gamma_total_pn_per_um",
                                     driving_correlation_time_s=DRIVING_TAU_S)
        assert report.verdict is StationarityVerdict.STATIONARY
        assert report.mean == pytest.approx(4.0, abs=0.02)
        naive = report.std / math.sqrt(report.n_eff * 2.0 * report.tau_int_s / dt)
        assert report.sem is not None and report.sem > naive, (
            "sem must exceed the naive sigma/sqrt(N); publishing the naive form is the standard way to "
            "quote an error bar that is far too small for a correlated series")
        assert report.n_eff < x.size

    def test_the_contract_is_recorded_next_to_the_verdict(self) -> None:
        """min_windows / drift_sigma are the caller's PRE-declared criterion and must survive into the record."""
        report = assess_stationarity(np.random.default_rng(7).normal(size=4_000), 0.01,
                                     observable="x", driving_correlation_time_s=DRIVING_TAU_S,
                                     min_windows=3.0, drift_sigma=2.0)
        assert report.contract["min_windows"] == 3.0
        assert report.contract["drift_sigma"] == 2.0
        assert report.contract["driving_correlation_time_s"] == pytest.approx(DRIVING_TAU_S)
        assert set(report.as_dict()) >= {"verdict", "mean", "sem", "tau_int_s", "n_eff", "contract"}

    def test_a_tighter_contract_can_only_refuse_more(self) -> None:
        """Monotonicity: raising min_windows must never turn a refusal into an acceptance.

        This is the property that makes the contract meaningful as a gate — if a stricter criterion could
        pass a series a looser one refused, declaring it in advance would protect nothing.
        """
        dt = 0.01
        x = 4.0 + _ar1(2_000, 0.9, seed=8, sigma=0.1)
        verdicts = [
            assess_stationarity(x, dt, observable="x", driving_correlation_time_s=DRIVING_TAU_S,
                                min_windows=mw).verdict
            for mw in (1.0, 5.0, 50.0)
        ]
        assert verdicts[0] is StationarityVerdict.STATIONARY
        assert verdicts[-1] is StationarityVerdict.TOO_SHORT

    def test_driving_time_must_be_positive_finite(self) -> None:
        """A zero or NaN driving time would make the required span zero — every series would pass."""
        x = np.random.default_rng(9).normal(size=64)
        with pytest.raises(ValueError, match="driving_correlation_time_s"):
            assess_stationarity(x, 0.01, observable="x", driving_correlation_time_s=0.0)

    def test_report_is_json_able(self) -> None:
        """``as_dict`` feeds an artifact measurements block, so every value must survive json.dumps."""
        import json

        report = assess_stationarity(np.random.default_rng(10).normal(size=4_000), 0.01,
                                     observable="x", driving_correlation_time_s=DRIVING_TAU_S)
        assert isinstance(report, StationarityReport)
        json.dumps(report.as_dict(), allow_nan=True)
