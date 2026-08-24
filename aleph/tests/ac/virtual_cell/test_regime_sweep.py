"""Analytic positive and negative controls for regime classification and the regime-aware sweep.

Every control here has a known answer by construction, so a misclassification is a measurement of the
classifier rather than an opinion about it.  The controls span the regimes the classifier claims to
separate plus three negative controls it must REFUSE:

============================  ==========================================================
control                       expected verdict
============================  ==========================================================
exponential relaxation        ``FIXED_POINT``
Ornstein-Uhlenbeck            ``STOCHASTIC_NESS``
supercritical Hopf, mu > 0    ``LIMIT_CYCLE``
subcritical Hopf, mu < 0      ``QUASICYCLE`` (noise-sustained)
two incommensurate tones      ``QUASIPERIODIC``
Lorenz, RK4-discretised       ``CHAOS``
logistic map, r = 3.99        ``CHAOS``
white noise                   ``UNDECIDED`` — no temporal structure to classify
linear drift + OU             ``UNDECIDED`` — a transient has no regime yet
300 samples                   ``UNDECIDED`` — below the declared minimum
============================  ==========================================================

ONE contract is declared for all of them, in :data:`CONTRACT`, before any of them is run.  It is not
tuned per control, and no threshold in it was chosen after seeing a verdict; where a threshold has an
empirical basis, that basis is a property of an estimator measured on the controls (the Hann window's
own linewidth, the Rosenstein estimator's noise floor) and is named in the comment beside it.

The suite also pins the two claims the module makes about its own provenance: that it agrees with
:mod:`aleph.engine.observe.stationarity` on the shared estimator, and that importing it does not
initialise Warp.
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
import textwrap
import warnings
from dataclasses import fields
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from aleph.virtual_cell.regime import (
    Regime,
    RegimeEvidenceContract,
    classify_regime,
    delay_embedding,
    embedding_delay,
    integrated_autocorrelation_time,
    normalized_autocorrelation,
    phase_randomized_surrogate,
    recurrence_determinism,
    spectral_peaks,
    welch_psd,
)
from aleph.virtual_cell.sweep import (
    EnsembleDefect,
    EnsembleMember,
    ExpansionRoute,
    SensitivityVerdict,
    SweepPoint,
    SweepVerdictKind,
    audit_ensemble,
    ess_interval,
    finite_difference_sensitivity,
    moment_closure_gap,
    phase_decomposition,
    summarize_ensemble,
)

# --------------------------------------------------------------------------------------------------
# The pre-declared evidence contract.  Declared once, before any control is run, and used for all.
# --------------------------------------------------------------------------------------------------
CONTRACT = RegimeEvidenceContract(
    # A Welch estimate with 4 segments needs a few hundred samples per segment to be worth reading.
    min_samples=512,
    # Matches the engine's stationarity default: 6 is the usual MD choice, 5 for shorter records.
    sokal_window_c=5.0,
    min_correlation_windows=20.0,
    min_n_eff=30.0,
    drift_sigma=1.0,
    fixed_point_relative_amplitude=0.01,
    welch_segment_count=4,
    # A chi-squared tail argument: with ~4 averaged segments a single bin exceeds 8x its continuum
    # with probability far below one per spectrum of ~600 bins.
    peak_prominence_ratio_min=8.0,
    baseline_octave_ratio=2.0,
    peak_power_fraction_min=0.05,
    secondary_peak_power_fraction=0.05,
    harmonic_relative_tolerance=0.02,
    max_rational_denominator=4,
    # A Hann window's own half-power main lobe is 1.44 resolution bins, so a resolution-limited line
    # cannot measure narrower than that; 2.5 allows a 70 % margin for discrete interpolation.
    limit_cycle_max_width_bins=2.5,
    # ~3x the window's own linewidth: broadened by the process, not by the window.
    quasicycle_min_width_bins=4.0,
    white_noise_tau_dt_max=1.0,
    white_noise_flatness_min=0.5,
    determinism_min=0.5,
    predictability_z_min=4.0,
    predictability_neighbours=4,
    n_surrogates=8,
    surrogate_seed=20260729,
    embedding_dimension=4,
    recurrence_rate_target=0.05,
    min_diagonal_length=2,
    recurrence_max_points=800,
    theiler_tau_multiple=2.0,
    lyapunov_fit_tau_multiple=1.0,
    lyapunov_min_fit_lags=5,
    lyapunov_max_fit_lags=60,
)

N_SAMPLES = 8000
DT = 0.05
COVERAGE_Z = 1.96


# --------------------------------------------------------------------------------------------------
# Analytic controls.
# --------------------------------------------------------------------------------------------------
def ornstein_uhlenbeck(
    n: int, dt: float, *, theta: float, sigma: float, seed: int, x0: float = 0.0
) -> np.ndarray:
    """Exact-update OU path: ``dx = -theta x dt + sigma dW``, stationary std ``sigma/sqrt(2 theta)``."""
    rng = np.random.default_rng(seed)
    decay = float(np.exp(-theta * dt))
    noise = sigma * float(np.sqrt((1.0 - decay * decay) / (2.0 * theta)))
    out = np.empty(n, dtype=np.float64)
    value = x0
    for i in range(n):
        value = decay * value + noise * rng.standard_normal()
        out[i] = value
    return out


def exponential_relaxation(n: int, dt: float, *, tau: float, seed: int, noise: float) -> np.ndarray:
    """A deterministic relaxation to a constant, plus a fluctuation far below the declared scale."""
    rng = np.random.default_rng(seed)
    t = np.arange(n, dtype=np.float64) * dt
    return 0.5 + 0.5 * np.exp(-t / tau) + noise * rng.standard_normal(n)


def hopf(n: int, dt: float, *, mu: float, omega: float, noise: float, seed: int) -> np.ndarray:
    """Supercritical Hopf normal form ``dz = ((mu + i w) z - |z|^2 z) dt + noise dW``; observe ``Re z``.

    ``mu > 0`` is above the bifurcation: a self-sustained limit cycle whose linewidth is set by the
    record length.  ``mu < 0`` is below it: a damped spiral kept alive by the noise, whose linewidth is
    set by ``2|mu|`` and is therefore genuinely broad.
    """
    rng = np.random.default_rng(seed)
    z = 0.1 + 0.0j
    out = np.empty(n, dtype=np.float64)
    root_dt = float(np.sqrt(dt))
    for i in range(n):
        drift = (mu + 1j * omega) * z - (abs(z) ** 2) * z
        z = z + drift * dt + noise * root_dt * (rng.standard_normal() + 1j * rng.standard_normal())
        out[i] = z.real
    return out


def two_tone(n: int, dt: float, *, f1: float, f2: float, noise: float, seed: int) -> np.ndarray:
    """Two incommensurate tones; ``f2/f1`` is the golden ratio, the least rationally approximable."""
    rng = np.random.default_rng(seed)
    t = np.arange(n, dtype=np.float64) * dt
    return (
        np.sin(2.0 * np.pi * f1 * t)
        + 0.8 * np.sin(2.0 * np.pi * f2 * t)
        + noise * rng.standard_normal(n)
    )


def lorenz_x(n: int, dt: float, *, transient: int = 2000) -> np.ndarray:
    """RK4-discretised Lorenz (sigma=10, rho=28, beta=8/3); observe ``x``.  True lambda_max ~0.9/s."""
    sigma, rho, beta = 10.0, 28.0, 8.0 / 3.0
    state = np.array([1.0, 1.0, 20.0], dtype=np.float64)

    def derivative(v: np.ndarray) -> np.ndarray:
        x, y, z = v
        return np.array([sigma * (y - x), x * (rho - z) - y, x * y - beta * z])

    out = np.empty(n, dtype=np.float64)
    for i in range(n + transient):
        k1 = derivative(state)
        k2 = derivative(state + 0.5 * dt * k1)
        k3 = derivative(state + 0.5 * dt * k2)
        k4 = derivative(state + dt * k3)
        state = state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        if i >= transient:
            out[i - transient] = state[0]
    return out


def logistic_map(n: int, *, r: float = 3.99, x0: float = 0.4, transient: int = 1000) -> np.ndarray:
    """Logistic map at ``r = 3.99``.  Spectrally white and delta-correlated: the ACF cannot see it."""
    value = x0
    out = np.empty(n, dtype=np.float64)
    for i in range(n + transient):
        value = r * value * (1.0 - value)
        if i >= transient:
            out[i - transient] = value
    return out


def white_noise(n: int, *, seed: int) -> np.ndarray:
    """Uncorrelated Gaussian samples — the negative control the classifier must refuse to label."""
    return np.random.default_rng(seed).standard_normal(n)


#: ``(name, series, dt, amplitude_scale, expected_regime)``.  Expectations fixed before the run.
CONTROLS: tuple[tuple[str, np.ndarray, float, float, Regime], ...] = (
    (
        "fixed_point",
        exponential_relaxation(N_SAMPLES, DT, tau=2.0, seed=11, noise=1e-4),
        DT,
        1.0,
        Regime.FIXED_POINT,
    ),
    (
        "ou_ness",
        ornstein_uhlenbeck(N_SAMPLES, DT, theta=0.8, sigma=0.5, seed=3),
        DT,
        1.0,
        Regime.STOCHASTIC_NESS,
    ),
    (
        "limit_cycle",
        hopf(N_SAMPLES, DT, mu=1.0, omega=2.0 * np.pi * 0.3, noise=0.02, seed=5),
        DT,
        1.0,
        Regime.LIMIT_CYCLE,
    ),
    (
        "quasicycle_weak",
        hopf(N_SAMPLES, DT, mu=-0.1, omega=2.0 * np.pi * 0.3, noise=0.06, seed=7),
        DT,
        1.0,
        Regime.QUASICYCLE,
    ),
    (
        "quasicycle_strong",
        hopf(N_SAMPLES, DT, mu=-0.25, omega=2.0 * np.pi * 0.3, noise=0.09, seed=7),
        DT,
        1.0,
        Regime.QUASICYCLE,
    ),
    (
        "quasiperiodic",
        two_tone(N_SAMPLES, DT, f1=0.11, f2=0.11 * 1.6180339887, noise=0.02, seed=9),
        DT,
        1.0,
        Regime.QUASIPERIODIC,
    ),
    ("lorenz", lorenz_x(N_SAMPLES, 0.02), 0.02, 10.0, Regime.CHAOS),
    ("logistic", logistic_map(N_SAMPLES), 1.0, 1.0, Regime.CHAOS),
    ("white_noise", white_noise(N_SAMPLES, seed=13), DT, 1.0, Regime.UNDECIDED),
    (
        "drifting_transient",
        np.arange(N_SAMPLES) * 0.001
        + ornstein_uhlenbeck(N_SAMPLES, DT, theta=0.8, sigma=0.2, seed=4),
        DT,
        1.0,
        Regime.UNDECIDED,
    ),
    (
        "too_short",
        ornstein_uhlenbeck(300, DT, theta=0.8, sigma=0.5, seed=3),
        DT,
        1.0,
        Regime.UNDECIDED,
    ),
)


@pytest.fixture(scope="module")
def classified() -> dict[str, Any]:
    """Classify every control once under the single declared contract."""
    results = {}
    for name, series, dt, scale, expected in CONTROLS:
        results[name] = {
            "report": classify_regime(
                series, dt, observable=name, amplitude_scale=scale, contract=CONTRACT
            ),
            "expected": expected,
            "series": series,
            "dt": dt,
        }
    return results


# --------------------------------------------------------------------------------------------------
# Regime classification.
# --------------------------------------------------------------------------------------------------
def test_control_confusion_matrix(classified: dict[str, Any], capsys: Any) -> None:
    """Every analytic control lands on its expected regime; the full matrix is printed."""
    labels = [
        Regime.FIXED_POINT,
        Regime.STOCHASTIC_NESS,
        Regime.QUASICYCLE,
        Regime.LIMIT_CYCLE,
        Regime.QUASIPERIODIC,
        Regime.CHAOS,
        Regime.UNDECIDED,
    ]
    matrix = {truth: dict.fromkeys(labels, 0) for truth in labels}
    rows = []
    for name, entry in classified.items():
        report = entry["report"]
        matrix[entry["expected"]][report.regime] += 1
        rows.append((name, entry["expected"].value, report.regime.value))

    with capsys.disabled():
        print("\n\nCONFUSION MATRIX (rows = truth, columns = predicted)")
        header = f"{'truth':>18s} | " + " ".join(f"{label.value[:9]:>9s}" for label in labels)
        print(header)
        print("-" * len(header))
        for truth in labels:
            counts = " ".join(f"{matrix[truth][pred]:>9d}" for pred in labels)
            print(f"{truth.value:>18s} | {counts}")
        print("\nper-control verdicts")
        for name, expected, actual in rows:
            flag = "ok " if expected == actual else "MISS"
            print(f"  {flag} {name:>20s}  expected {expected:<16s} got {actual}")
        print("\nkey statistics")
        for name, entry in classified.items():
            stats = entry["report"].statistics
            print(
                f"  {name:>20s} tau/dt={stats.tau_int_s / stats.dt_s:8.3g} "
                f"n_eff={stats.n_eff:9.4g} flat={stats.spectral_flatness:9.3g} "
                f"fund={len(stats.fundamental_frequencies_hz)} det={stats.determinism:6.3g} "
                f"pred_z={stats.predictability_z:7.3g} lam={stats.lyapunov_per_s:8.3g}/s"
            )

    wrong = [(n, e, a) for n, e, a in rows if e != a]
    assert not wrong, f"misclassified controls: {wrong}"
    for truth in labels:
        assert matrix[truth][truth] == sum(matrix[truth].values())


def test_undecided_is_a_real_and_frequent_answer(classified: dict[str, Any]) -> None:
    """Three distinct negative controls refuse, each for its own named reason."""
    undecided = {
        name: entry["report"]
        for name, entry in classified.items()
        if entry["report"].regime is Regime.UNDECIDED
    }
    assert set(undecided) == {"white_noise", "drifting_transient", "too_short"}
    reasons = {name: report.reason for name, report in undecided.items()}
    assert "indistinguishable" in reasons["white_noise"]
    assert "transient" in reasons["drifting_transient"]
    assert "minimum" in reasons["too_short"]
    for report in undecided.values():
        assert report.expansions, "an UNDECIDED verdict must name where it routes"
        assert not report.is_decided


@pytest.mark.parametrize("magnitude", [1e180, 1e200, 1e250])
def test_an_overflowing_second_moment_is_undecided_and_not_a_traceback(magnitude: float) -> None:
    """Finite samples whose SECOND moment overflows must yield UNDECIDED, not a ValueError.

    Regression for the crash the S0 invariance battery hit.  ``_as_series`` guarantees every
    SAMPLE is finite, which is not the same as every STATISTIC being finite: an explicit scheme
    one step below overflow produces samples around 1e200 that are all valid doubles, and then
    the variance overflows to ``inf`` and ``tau_int`` becomes ``nan``.  That used to reach
    ``round(nan)`` in the Theiler window and raise ``cannot convert float NaN to integer``, sixty
    lines from the cause and reading like a library bug.  Handing a classifier a diverged run is
    legitimate — the divergence IS the finding — so the reply must be a verdict.

    The input guard cannot catch this: every sample here passes ``np.isfinite``.
    """
    rng = np.random.default_rng(3)
    series = rng.standard_normal(3000) * magnitude
    assert np.isfinite(series).all(), "the point of this test is that the SAMPLES are finite"

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        report = classify_regime(
            series, DT, observable="overflowed", amplitude_scale=1.0, contract=CONTRACT
        )
    assert report.regime is Regime.UNDECIDED
    assert not report.is_decided
    assert "not finite" in report.reason
    assert report.expansions, "a divergence must route somewhere"


def test_a_well_scaled_series_is_still_classified_after_the_overflow_guard() -> None:
    """The guard must not swallow ordinary series — a large but non-overflowing scale still decides."""
    series = ornstein_uhlenbeck(4000, DT, theta=0.8, sigma=0.5, seed=11) * 1e6
    report = classify_regime(
        series, DT, observable="large_but_finite", amplitude_scale=1e6, contract=CONTRACT
    )
    assert report.regime is not Regime.UNDECIDED or "not finite" not in report.reason


def test_nonzero_fluctuation_is_never_scored_as_failure(classified: dict[str, Any]) -> None:
    """A NESS with a large sustained fluctuation is an accepted state, not a convergence failure."""
    report = classified["ou_ness"]["report"]
    assert report.regime is Regime.STOCHASTIC_NESS
    assert report.is_decided
    assert report.statistics.std > 0.3
    assert "not a convergence failure" in report.reason


def test_oscillatory_regimes_flag_that_a_mean_is_incomplete(classified: dict[str, Any]) -> None:
    """The three orbit regimes report ``is_oscillatory``; the two non-orbit ones do not."""
    oscillatory = {"limit_cycle", "quasicycle_weak", "quasicycle_strong", "quasiperiodic"}
    for name, entry in classified.items():
        assert entry["report"].is_oscillatory is (name in oscillatory), name
    period = classified["limit_cycle"]["report"].dominant_period_s
    assert period is not None
    assert period == pytest.approx(1.0 / 0.3, rel=0.05)


def test_regime_report_is_an_unverified_observer(classified: dict[str, Any]) -> None:
    """No report may claim any evidence authority beyond ``unverified-observer``."""
    for entry in classified.values():
        payload = entry["report"].as_dict()
        assert payload["evidence_authority"] == "unverified-observer"
        assert "native" not in payload["evidence_authority"]


def test_contract_has_no_defaults_and_forces_an_undecided_band() -> None:
    """Every threshold is the caller's, and no contract may partition sharp against broad."""
    import dataclasses

    for descriptor in fields(RegimeEvidenceContract):
        assert descriptor.default is dataclasses.MISSING, descriptor.name
        assert descriptor.default_factory is dataclasses.MISSING, descriptor.name

    payload = {
        descriptor.name: getattr(CONTRACT, descriptor.name)
        for descriptor in fields(RegimeEvidenceContract)
    }
    payload["quasicycle_min_width_bins"] = payload["limit_cycle_max_width_bins"]
    with pytest.raises(ValueError, match="UNDECIDED band"):
        RegimeEvidenceContract(**payload)


def test_classifier_rejects_an_undeclared_contract() -> None:
    """A caller cannot pass a bare mapping and have thresholds invented for it."""
    with pytest.raises(TypeError, match="declared before the run"):
        classify_regime(
            np.zeros(1024),
            DT,
            observable="x",
            amplitude_scale=1.0,
            contract={"min_samples": 10},  # type: ignore[arg-type]
        )


def test_theiler_exclusion_removes_trivial_determinism() -> None:
    """Without the Theiler window a smooth stochastic path is spuriously 'deterministic'."""
    series = ornstein_uhlenbeck(4000, DT, theta=0.8, sigma=0.5, seed=3)
    acf = normalized_autocorrelation(series)
    embedded = delay_embedding(series, delay=embedding_delay(acf, max_delay=500), dimension=4)
    _, without, _ = recurrence_determinism(
        embedded,
        recurrence_rate_target=0.05,
        min_diagonal_length=2,
        theiler_window=1,
        max_points=800,
    )
    _, with_window, _ = recurrence_determinism(
        embedded,
        recurrence_rate_target=0.05,
        min_diagonal_length=2,
        theiler_window=200,
        max_points=800,
    )
    assert without > with_window
    assert 0.0 <= with_window <= 1.0


def test_surrogate_preserves_the_power_spectrum() -> None:
    """The null must be a spectrum-matched null, or the determinism test is a spectrum test."""
    series = lorenz_x(2000, 0.02)
    surrogate = phase_randomized_surrogate(series, np.random.default_rng(0))
    original = np.abs(np.fft.rfft(series - series.mean()))
    matched = np.abs(np.fft.rfft(surrogate - surrogate.mean()))
    np.testing.assert_allclose(original, matched, rtol=1e-9, atol=1e-8)
    assert not np.allclose(series, surrogate)


def test_peak_width_is_normalized_by_resolution_not_by_record_length() -> None:
    """A pure tone measures the same width in BINS at two record lengths, though its Hz width halves."""
    widths, hz = [], []
    for n in (4000, 8000):
        t = np.arange(n) * DT
        series = np.sin(2.0 * np.pi * 0.3 * t)
        freq, psd, resolution = welch_psd(series, DT, n_segments=4)
        peaks = spectral_peaks(
            freq,
            psd,
            resolution,
            prominence_ratio_min=8.0,
            baseline_octave_ratio=2.0,
            power_fraction_min=0.05,
        )
        assert peaks
        widths.append(peaks[0].width_bins)
        hz.append(peaks[0].width_hz)
    assert widths[1] == pytest.approx(widths[0], rel=0.3)
    assert hz[1] < 0.75 * hz[0]


def test_integrated_autocorrelation_time_matches_the_engine_estimator() -> None:
    """The re-expressed estimator must equal the engine's, not merely resemble it.

    Loaded by file path so the comparison does not import ``aleph.engine``, whose package import
    initialises Warp — which this CPU-only layer forbids.
    """
    path = Path(__file__).resolve().parents[4] / "aleph/engine/observe/stationarity.py"
    if not path.exists():  # pragma: no cover - the engine module is expected to be present
        pytest.skip(f"engine stationarity module not found at {path}")
    spec = importlib.util.spec_from_file_location("_engine_stationarity_oracle", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
        series = ornstein_uhlenbeck(3000, DT, theta=0.8, sigma=0.5, seed=17)
        engine_tau, engine_n_eff, _ = module.integrated_autocorrelation_time(series, DT, c=5.0)
        mine_tau, mine_n_eff, _ = integrated_autocorrelation_time(series, DT, sokal_window_c=5.0)
        assert mine_tau == pytest.approx(engine_tau, rel=1e-12)
        assert mine_n_eff == pytest.approx(engine_n_eff, rel=1e-12)
    finally:
        sys.modules.pop(spec.name, None)


def test_importing_this_layer_does_not_initialise_warp() -> None:
    """The CPU-only contract, asserted rather than assumed — in a SUBPROCESS.

    Asserting against this process's ``sys.modules`` measures whether ANY test in the whole
    pytest session ever imported Warp, which is a different claim from the one the name makes
    and a false one: ``tests/ac`` contains tests that legitimately import it.  The in-process
    form therefore passed when this file ran alone and failed under the full suite.  A fresh
    interpreter measures what the name says, and it covers the file-path engine load too — the
    claim being that reaching the engine's estimator does not drag the engine package in.
    """
    script = textwrap.dedent(
        """
        import importlib.util
        import sys
        from pathlib import Path

        from aleph.virtual_cell.regime import integrated_autocorrelation_time
        from aleph.virtual_cell.sweep import finite_difference_sensitivity  # noqa: F401

        path = Path("aleph/engine/observe/stationarity.py").resolve()
        spec = importlib.util.spec_from_file_location("_engine_stationarity_oracle", path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)

        assert "warp" not in sys.modules, "the CPU-only layer initialised Warp"
        assert "aleph.engine" not in sys.modules, "file-path load dragged in the package"
        print("clean")
        """
    )
    done = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        cwd=str(Path(__file__).resolve().parents[4]),
    )
    assert done.returncode == 0, done.stderr
    assert "clean" in done.stdout


# --------------------------------------------------------------------------------------------------
# ESS-aware error bars.
# --------------------------------------------------------------------------------------------------
def autoregressive(n: int, *, phi: float, seed: int) -> np.ndarray:
    """AR(1) with a known correlation time ``tau_int/dt = (1 + phi) / (2 (1 - phi))``."""
    rng = np.random.default_rng(seed)
    out = np.empty(n, dtype=np.float64)
    value = 0.0
    scale = float(np.sqrt(1.0 - phi * phi))
    for i in range(n):
        value = phi * value + scale * rng.standard_normal()
        out[i] = value
    return out


def test_strong_autocorrelation_collapses_ess_and_widens_the_interval(capsys: Any) -> None:
    """``n_eff`` must fall far below ``N``, and the widened interval is what gets reported."""
    series = autoregressive(8000, phi=0.98, seed=101)
    interval = ess_interval(
        series, DT, observable="correlated", sokal_window_c=5.0, coverage_z=COVERAGE_Z
    )
    expected_tau_samples = (1.0 + 0.98) / (2.0 * (1.0 - 0.98))
    with capsys.disabled():
        print(
            f"\nESS: N={interval.n_samples} n_eff={interval.n_eff:.4g} "
            f"tau/dt={interval.tau_int_s / DT:.4g} (analytic {expected_tau_samples:.4g}) "
            f"sem_naive={interval.sem_naive:.5g} sem_ess={interval.sem_ess:.5g} "
            f"widening={interval.widening_factor:.4g}x"
        )
    assert interval.n_eff < 0.05 * interval.n_samples
    assert interval.tau_int_s / DT == pytest.approx(expected_tau_samples, rel=0.35)
    assert interval.widening_factor > 5.0
    assert interval.high - interval.low == pytest.approx(2.0 * COVERAGE_Z * interval.sem_ess)
    assert interval.sem_ess > interval.sem_naive


def test_uncorrelated_series_keeps_the_naive_width() -> None:
    """The ESS form must not be a blanket inflation; on white noise it collapses onto ``sigma/sqrt(N)``."""
    interval = ess_interval(
        white_noise(8000, seed=5), DT, observable="iid", sokal_window_c=5.0, coverage_z=1.0
    )
    assert interval.widening_factor == pytest.approx(1.0, abs=0.05)
    assert interval.n_eff == pytest.approx(interval.n_samples, rel=0.05)


# --------------------------------------------------------------------------------------------------
# Sweep sensitivity: the cross-regime refusal.
# --------------------------------------------------------------------------------------------------
def make_point(series: np.ndarray, dt: float, *, value: float, scale: float = 1.0) -> SweepPoint:
    """Build a sweep point at parameter ``value`` under the module contract."""
    return SweepPoint.measure(
        series,
        dt,
        parameter_name="mu",
        parameter_value=value,
        parameter_unit="1",
        observable="edge_position",
        observable_unit="um",
        amplitude_scale=scale,
        regime_contract=CONTRACT,
        coverage_z=COVERAGE_Z,
    )


@pytest.fixture(scope="module")
def bifurcation_pair() -> tuple[SweepPoint, SweepPoint]:
    """Two points either side of a Hopf boundary: a NESS below, a limit cycle above."""
    below = make_point(
        ornstein_uhlenbeck(N_SAMPLES, DT, theta=0.8, sigma=0.5, seed=3), DT, value=-0.5
    )
    above = make_point(
        hopf(N_SAMPLES, DT, mu=1.0, omega=2.0 * np.pi * 0.3, noise=0.02, seed=5), DT, value=1.0
    )
    return below, above


def test_sensitivity_is_refused_across_a_regime_boundary(
    bifurcation_pair: tuple[SweepPoint, SweepPoint], capsys: Any
) -> None:
    """Two points in different regimes must not yield a number."""
    below, above = bifurcation_pair
    assert below.regime.regime is Regime.STOCHASTIC_NESS
    assert above.regime.regime is Regime.LIMIT_CYCLE
    verdict = finite_difference_sensitivity(below, above, min_n_eff=30.0)
    with capsys.disabled():
        print(f"\ncross-regime verdict: {verdict.kind.value}")
        print(f"  reason: {verdict.reason}")
        print(f"  expansions: {[route.value for route in verdict.expansions]}")
    assert verdict.kind is SweepVerdictKind.REFUSED_CROSS_REGIME
    assert verdict.value is None
    assert verdict.uncertainty is None
    assert not verdict.is_reported
    assert ExpansionRoute.REGIME_PARTITION in verdict.expansions
    assert ExpansionRoute.LATENT_COORDINATE in verdict.expansions


def test_sensitivity_is_reported_within_one_regime(capsys: Any) -> None:
    """Two NESS points at different parameter values do yield a number, with an ESS uncertainty."""
    low = make_point(ornstein_uhlenbeck(N_SAMPLES, DT, theta=0.8, sigma=0.5, seed=3), DT, value=1.0)
    high = make_point(
        ornstein_uhlenbeck(N_SAMPLES, DT, theta=0.8, sigma=0.5, seed=23) + 0.4, DT, value=2.0
    )
    verdict = finite_difference_sensitivity(low, high, min_n_eff=30.0)
    with capsys.disabled():
        print(
            f"\nwithin-regime sensitivity: {verdict.value:.5g} +/- {verdict.uncertainty:.5g} "
            f"{verdict.unit} (z={verdict.coverage_z})"
        )
    assert verdict.kind is SweepVerdictKind.REPORTED
    assert verdict.is_reported
    assert verdict.unit == "um/1"
    assert verdict.value == pytest.approx(0.4, abs=0.2)
    naive = np.hypot(low.interval.sem_naive, high.interval.sem_naive) * COVERAGE_Z
    assert verdict.uncertainty > naive


def test_sensitivity_is_refused_when_a_point_is_undecided() -> None:
    """An undecided regime cannot be differenced against anything."""
    decided = make_point(
        ornstein_uhlenbeck(N_SAMPLES, DT, theta=0.8, sigma=0.5, seed=3), DT, value=1.0
    )
    unlabelled = make_point(white_noise(N_SAMPLES, seed=13), DT, value=2.0)
    assert unlabelled.regime.regime is Regime.UNDECIDED
    verdict = finite_difference_sensitivity(decided, unlabelled, min_n_eff=30.0)
    assert verdict.kind is SweepVerdictKind.REFUSED_UNDECIDED_REGIME
    assert verdict.value is None
    assert ExpansionRoute.TIME_MODALITY in verdict.expansions


def test_sensitivity_is_refused_when_the_effective_sample_size_is_too_small() -> None:
    """The step count is not the sample size; the floor is checked against ``n_eff``."""
    low = make_point(ornstein_uhlenbeck(N_SAMPLES, DT, theta=0.8, sigma=0.5, seed=3), DT, value=1.0)
    high = make_point(
        ornstein_uhlenbeck(N_SAMPLES, DT, theta=0.8, sigma=0.5, seed=23), DT, value=2.0
    )
    assert low.interval.n_samples == N_SAMPLES
    assert low.interval.n_eff < 0.1 * N_SAMPLES
    verdict = finite_difference_sensitivity(low, high, min_n_eff=1000.0)
    assert verdict.kind is SweepVerdictKind.REFUSED_INSUFFICIENT_ESS
    assert verdict.value is None
    assert ExpansionRoute.TIME_MODALITY in verdict.expansions


def test_sensitivity_rejects_a_degenerate_or_mismatched_pair(
    bifurcation_pair: tuple[SweepPoint, SweepPoint],
) -> None:
    """A zero separation or a mismatched observable is a caller error, not a finding."""
    below, _ = bifurcation_pair
    with pytest.raises(ValueError, match="nonzero"):
        finite_difference_sensitivity(below, below, min_n_eff=30.0)
    other = SweepPoint.measure(
        below.series,
        DT,
        parameter_name="mu",
        parameter_value=2.0,
        parameter_unit="1",
        observable="tension",
        observable_unit="pN/um",
        amplitude_scale=1.0,
        regime_contract=CONTRACT,
        coverage_z=COVERAGE_Z,
    )
    with pytest.raises(ValueError, match="same observable"):
        finite_difference_sensitivity(below, other, min_n_eff=30.0)


def test_a_refusal_must_name_an_expansion() -> None:
    """A refusal with no route is how a hard problem gets dropped instead of scheduled."""
    with pytest.raises(ValueError, match="expansion route"):
        SensitivityVerdict(
            kind=SweepVerdictKind.REFUSED_CROSS_REGIME,
            parameter_name="mu",
            observable="x",
            unit="um/1",
            low_value=0.0,
            high_value=1.0,
            low_regime=Regime.STOCHASTIC_NESS,
            high_regime=Regime.LIMIT_CYCLE,
            reason="no route given",
        )
    with pytest.raises(ValueError, match="must not carry a value"):
        SensitivityVerdict(
            kind=SweepVerdictKind.REFUSED_CROSS_REGIME,
            parameter_name="mu",
            observable="x",
            unit="um/1",
            low_value=0.0,
            high_value=1.0,
            low_regime=Regime.STOCHASTIC_NESS,
            high_regime=Regime.LIMIT_CYCLE,
            value=1.0,
            expansions=(ExpansionRoute.REGIME_PARTITION,),
        )


# --------------------------------------------------------------------------------------------------
# Phase / period decomposition.
# --------------------------------------------------------------------------------------------------
def test_phase_decomposition_recovers_what_the_time_average_deleted(capsys: Any) -> None:
    """A limit cycle's mean is near zero while its phase profile spans the full orbit."""
    series = hopf(N_SAMPLES, DT, mu=1.0, omega=2.0 * np.pi * 0.3, noise=0.02, seed=5)
    report = classify_regime(
        series, DT, observable="edge_position", amplitude_scale=1.0, contract=CONTRACT
    )
    assert report.regime is Regime.LIMIT_CYCLE
    period = report.dominant_period_s
    assert period is not None
    decomposition = phase_decomposition(
        series, DT, observable="edge_position", period_s=period, n_bins=16
    )
    with capsys.disabled():
        print(
            f"\nphase decomposition: period={decomposition.period_s:.4g} s over "
            f"{decomposition.n_cycles:.3g} cycles; time_average={decomposition.time_average:.4g}, "
            f"cycle_amplitude={decomposition.cycle_amplitude:.4g} "
            f"({decomposition.amplitude_over_sem:.3g} sigma)"
        )
    assert abs(decomposition.time_average) < 0.1 * decomposition.cycle_amplitude
    assert decomposition.cycle_amplitude > 1.0
    assert decomposition.amplitude_over_sem > 3.0
    assert int(decomposition.phase_bin_counts.sum()) == series.size
    assert decomposition.phase_bin_means.shape == (16,)


def test_phase_decomposition_refuses_less_than_one_cycle() -> None:
    """A fold over a partial cycle is not a phase average."""
    with pytest.raises(ValueError, match="less than one"):
        phase_decomposition(np.zeros(100), DT, observable="x", period_s=100.0, n_bins=8)
    with pytest.raises(ValueError, match="cannot fill"):
        phase_decomposition(np.zeros(6), DT, observable="x", period_s=0.05, n_bins=8)


# --------------------------------------------------------------------------------------------------
# The three ensemble defects.
# --------------------------------------------------------------------------------------------------
def clean_members() -> list[EnsembleMember]:
    """Four well-formed, independently seeded members."""
    return [
        EnsembleMember(
            member_id=f"m{i}",
            observables={"tension": 3.0 + 0.1 * i, "area": 900.0 + i},
            seed=100 + i,
        )
        for i in range(4)
    ]


def test_missing_observable_is_detected_not_averaged_over() -> None:
    """Defect (a): a member lacking an observable its siblings have."""
    members = clean_members()
    members[2] = EnsembleMember(member_id="m2", observables={"area": 902.0}, seed=102)
    defects = audit_ensemble(members, weight_tolerance=1e-9)
    kinds = {record.defect for record in defects}
    assert EnsembleDefect.MISSING_OBSERVABLE in kinds
    record = next(r for r in defects if r.defect is EnsembleDefect.MISSING_OBSERVABLE)
    assert record.member_ids == ("m2",)
    summary = summarize_ensemble(members, weight_tolerance=1e-9)
    assert not summary.is_clean
    assert summary.means == {}, "a defective ensemble must not emit a partial mean"


def test_duplicate_seed_is_detected_and_deflates_the_independent_count() -> None:
    """Defect (b): two members that ran the same stream are not independent draws."""
    members = clean_members()
    members[3] = EnsembleMember(
        member_id="m3", observables={"tension": 3.3, "area": 903.0}, seed=100
    )
    defects = audit_ensemble(members, weight_tolerance=1e-9)
    record = next(r for r in defects if r.defect is EnsembleDefect.DUPLICATE_SEED)
    assert set(record.member_ids) == {"m0", "m3"}
    summary = summarize_ensemble(members, weight_tolerance=1e-9)
    assert summary.apparent_n == 4
    assert summary.independent_n == 3
    assert not summary.is_clean


def test_missing_seed_is_its_own_defect() -> None:
    """An unrecorded seed is unknown independence, which is not the same as independence."""
    members = clean_members()
    members[1] = EnsembleMember(member_id="m1", observables={"tension": 3.1, "area": 901.0})
    defects = audit_ensemble(members, weight_tolerance=1e-9)
    record = next(r for r in defects if r.defect is EnsembleDefect.MISSING_SEED)
    assert record.member_ids == ("m1",)


def test_mixture_weights_must_be_nonnegative_and_normalized() -> None:
    """Defect (c): weights that do not sum to one, or contain a negative, are rejected."""
    base = clean_members()
    unnormalized = [
        EnsembleMember(m.member_id, m.observables, m.seed, weight=w)
        for m, w in zip(base, (0.3, 0.3, 0.3, 0.3), strict=True)
    ]
    kinds = {r.defect for r in audit_ensemble(unnormalized, weight_tolerance=1e-9)}
    assert EnsembleDefect.WEIGHTS_NOT_NORMALIZED in kinds

    negative = [
        EnsembleMember(m.member_id, m.observables, m.seed, weight=w)
        for m, w in zip(base, (0.6, 0.6, -0.1, -0.1), strict=True)
    ]
    kinds = {r.defect for r in audit_ensemble(negative, weight_tolerance=1e-9)}
    assert EnsembleDefect.NEGATIVE_WEIGHT in kinds

    partial = [
        EnsembleMember(base[0].member_id, base[0].observables, base[0].seed, weight=0.5),
        *base[1:],
    ]
    kinds = {r.defect for r in audit_ensemble(partial, weight_tolerance=1e-9)}
    assert EnsembleDefect.PARTIAL_WEIGHTS in kinds


def test_a_clean_ensemble_summarizes_unweighted_and_weighted() -> None:
    """The audit is not a blanket refusal: a well-formed ensemble averages, weights included."""
    members = clean_members()
    summary = summarize_ensemble(members, weight_tolerance=1e-9)
    assert summary.is_clean
    assert summary.independent_n == 4
    assert summary.means["tension"] == pytest.approx(3.15)
    assert not summary.weighted

    weighted = [
        EnsembleMember(m.member_id, m.observables, m.seed, weight=w)
        for m, w in zip(members, (0.4, 0.3, 0.2, 0.1), strict=True)
    ]
    summary = summarize_ensemble(weighted, weight_tolerance=1e-9)
    assert summary.is_clean and summary.weighted
    expected = 0.4 * 3.0 + 0.3 * 3.1 + 0.2 * 3.2 + 0.1 * 3.3
    assert summary.means["tension"] == pytest.approx(expected)


def test_audit_reports_every_defect_at_once() -> None:
    """A caller handed the whole list fixes the ensemble once, not once per rerun."""
    members = [
        EnsembleMember("m0", {"tension": 3.0, "area": 900.0}, seed=7, weight=0.5),
        EnsembleMember("m1", {"area": 901.0}, seed=7, weight=0.9),
        EnsembleMember("m2", {"tension": 3.2, "area": 902.0}, weight=-0.1),
    ]
    kinds = {record.defect for record in audit_ensemble(members, weight_tolerance=1e-9)}
    assert kinds == {
        EnsembleDefect.MISSING_OBSERVABLE,
        EnsembleDefect.DUPLICATE_SEED,
        EnsembleDefect.MISSING_SEED,
        EnsembleDefect.NEGATIVE_WEIGHT,
        EnsembleDefect.WEIGHTS_NOT_NORMALIZED,
    }


# --------------------------------------------------------------------------------------------------
# The Jensen term as a measured moment-closure cost.
# --------------------------------------------------------------------------------------------------
def test_closure_gap_for_the_square_equals_the_variance_exactly(capsys: Any) -> None:
    """``E[X^2] - (E[X])^2 = Var[X]`` is an identity, so this is a check against truth."""
    samples = np.random.default_rng(2026).normal(loc=4.0, scale=1.5, size=200000)
    report = moment_closure_gap(samples, lambda x: x**2, observable="x_squared", tolerance=1.0)
    population_variance = float(np.var(samples))
    with capsys.disabled():
        print(
            f"\nclosure gap (f=x^2): E[f(X)]={report.ensemble_mean_of_f:.6f} "
            f"f(E[X])={report.f_of_ensemble_mean:.6f} gap={report.gap:.6f} "
            f"population Var={population_variance:.6f} (analytic sigma^2=2.25) "
            f"relative_gap={report.relative_gap:.5g}"
        )
    assert report.gap == pytest.approx(population_variance, rel=1e-12)
    assert report.gap == pytest.approx(2.25, rel=0.02)
    assert report.within_tolerance


def test_closure_gap_for_the_exponential_matches_the_lognormal_closed_form(capsys: Any) -> None:
    """``E[e^X] - e^{E[X]} = e^{mu}(e^{sigma^2/2} - 1)`` for Gaussian ``X``."""
    mu, sigma = 0.5, 0.4
    samples = np.random.default_rng(77).normal(loc=mu, scale=sigma, size=400000)
    report = moment_closure_gap(samples, np.exp, observable="exp_x", tolerance=1.0)
    analytic = float(np.exp(mu) * (np.exp(sigma**2 / 2.0) - 1.0))
    with capsys.disabled():
        print(
            f"closure gap (f=exp): measured={report.gap:.6f} analytic={analytic:.6f} "
            f"relative error={abs(report.gap - analytic) / analytic:.3e}"
        )
    assert report.gap == pytest.approx(analytic, rel=0.01)


def test_closure_gap_vanishes_for_a_linear_observable() -> None:
    """The gap is a property of nonlinearity; a linear observable must measure exactly zero."""
    samples = np.random.default_rng(5).normal(size=50000)
    report = moment_closure_gap(
        samples, lambda x: 3.0 * x - 2.0, observable="linear", tolerance=1e-12
    )
    assert report.gap == pytest.approx(0.0, abs=1e-10)
    assert report.within_tolerance
    assert report.expansions == ()


def test_closure_gap_breach_routes_to_a_richer_representation() -> None:
    """A gap past tolerance expands the representation; it never indicts the reference dynamics."""
    samples = np.random.default_rng(9).normal(loc=0.0, scale=3.0, size=20000)
    report = moment_closure_gap(samples, lambda x: x**2, observable="x_squared", tolerance=0.01)
    assert not report.within_tolerance
    assert ExpansionRoute.ENSEMBLE in report.expansions
    payload = report.as_dict()
    assert "NOT a defect in the reference dynamics" in payload["interpretation"]


def test_closure_gap_honours_and_validates_mixture_weights() -> None:
    """A weighted measure is supported; a malformed one is refused rather than renormalised."""
    samples = np.array([0.0, 2.0, 4.0])
    weights = np.array([0.25, 0.5, 0.25])
    report = moment_closure_gap(
        samples, lambda x: x**2, observable="x_squared", tolerance=1.0, weights=weights
    )
    mean = float(np.dot(weights, samples))
    expected = float(np.dot(weights, samples**2)) - mean**2
    assert report.gap == pytest.approx(expected, rel=1e-12)
    assert report.weighted
    with pytest.raises(ValueError, match="sum to 1"):
        moment_closure_gap(
            samples,
            lambda x: x**2,
            observable="q",
            tolerance=1.0,
            weights=np.array([0.5, 0.5, 0.5]),
        )
    with pytest.raises(ValueError, match="nonnegative"):
        moment_closure_gap(
            samples,
            lambda x: x**2,
            observable="q",
            tolerance=1.0,
            weights=np.array([1.2, -0.1, -0.1]),
        )
