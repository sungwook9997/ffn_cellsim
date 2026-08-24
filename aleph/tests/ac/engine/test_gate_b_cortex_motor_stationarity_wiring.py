"""The GATE-B cortex-motor driver must judge its own settling, and know its own bound-head lifetime.

WHY.  On 2026-07-28 a six-point blebbistatin sweep ran this driver for a FIXED 60 steps per point and
came out non-monotonic in ``k_xb``.  That was a protocol error, not physics: a softer crossbridge
relaxes more slowly, so a fixed step count compares the points at four different places on four
different transients.  ``observe/stationarity.py`` was written in the same session to make that a
verdict and was not wired into the driver it exists for.  These tests pin the wiring.

The second half is the driving time itself.  The driver stamped ``bound_head_lifetime_s = 1/0.4`` as a
literal, and every number in it was wrong for the runs it stamped: this repo's ``NMII_KOFF0`` is 0.35,
the sweep ran ``--catch-slip`` whose zero-load off-rate is ``k_catch0 + k_slip0`` and not the Bell
prefactor at all, and heads under load on a CATCH bond live longer than either.  That constant is the
denominator the whole contract is expressed in, so it is tested rather than trusted.

The helpers moved OUT of the driver on 2026-07-28 — into
:func:`aleph.components.motor.bell_kinetics_analytic.bound_state_lifetime` and
:func:`aleph.engine.observe.stationarity.assess_observables` — because a driver-local copy is the
defect this file documents: it is how the wrong literal survived in one driver while
``ac_gate_b_sf_motor_native.py`` could not express the question at all (it mentions a bound-head lifetime
zero times).  These tests exercise the driver's thin wrappers, which is what a caller actually touches,
and one case pins the shared function directly so a second lane inherits the fixed version.

CPU-only: both helpers are host-side analysis over arrays a run already produced.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from aleph.engine.cortex_motor_slice import CortexMotorParams
from aleph.engine.observe.stationarity import StationarityVerdict
from aleph.components.motor.segment_motor import SegmentDetachKinetics
from aleph.components.motor.bell_kinetics_analytic import bound_state_lifetime
from aleph.scripts.ac_gate_b_cortex_motor_native import (
    CATCH_SLIP_PROVISIONAL,
    STATIONARITY_OBSERVABLES,
    _assess_series,
    _bound_head_lifetime,
)

_BASE = dict(v0=0.5, f_stall=0.5, kappa=0.5, k_xb=1000.0, r0_head=0.02, r0_xb=0.0,
             capture_radius=0.05, k_on=50.0, k_off0=0.35, f0=7.13)


class _FakeArray:
    """Minimal stand-in for a Warp array: the helpers only ever call ``.numpy()`` on these."""

    def __init__(self, values: np.ndarray) -> None:
        self._values = values

    def numpy(self) -> np.ndarray:
        return self._values


class _FakeSlice:
    """A slice whose connector state carries a chosen binding pattern and per-head Bell load."""

    def __init__(self, bound: np.ndarray, loads: np.ndarray) -> None:
        state = type("S", (), {"bound_d": _FakeArray(bound), "loads_bell_d": _FakeArray(loads)})()
        self.connector = type("C", (), {"state": state})()


def _slip_params() -> CortexMotorParams:
    return CortexMotorParams(detach_kinetics=SegmentDetachKinetics.SLIP, **_BASE)


def _catch_slip_params() -> CortexMotorParams:
    return CortexMotorParams(detach_kinetics=SegmentDetachKinetics.CATCH_SLIP,
                             **CATCH_SLIP_PROVISIONAL, **_BASE)


def _trajectory(gamma: np.ndarray) -> list[dict]:
    """Wrap a gamma series in the row shape the driver records, filling the other observables."""
    return [{"gamma_total_pn_per_um": float(g), "gamma_network_pn_per_um": float(g) * 0.99,
             "gamma_source_pn_per_um": float(g) * 0.01, "bound_fraction": 0.98,
             "step": i, "t_s": (i + 1) * 0.01} for i, g in enumerate(gamma)]


class TestBoundHeadLifetime:
    """The driving correlation time is derived from the kinetics the run actually used."""

    def test_slip_zero_load_is_the_bell_prefactor(self) -> None:
        """In SLIP mode the zero-load lifetime is ``1/k_off0`` — 1/0.35 s, not the stamped 1/0.4."""
        result = _bound_head_lifetime(_FakeSlice(np.zeros(4, np.int32), np.zeros(4)), _slip_params())
        assert result["zero_load_s"] == pytest.approx(1.0 / 0.35)
        assert result["zero_load_s"] != pytest.approx(1.0 / 0.4, rel=1e-3)

    def test_catch_slip_zero_load_is_the_SUM_of_the_two_pathways(self) -> None:
        """Pereverzev at zero load is ``k_catch0 + k_slip0``; using either alone doubles the lifetime.

        This is the constant the 2026-07-28 sweep was judged against, and the sweep ran ``--catch-slip``.
        """
        result = _bound_head_lifetime(_FakeSlice(np.zeros(4, np.int32), np.zeros(4)),
                                      _catch_slip_params())
        expected = 1.0 / (CATCH_SLIP_PROVISIONAL["k_catch0"] + CATCH_SLIP_PROVISIONAL["k_slip0"])
        assert result["zero_load_s"] == pytest.approx(expected)
        assert result["zero_load_s"] == pytest.approx(1.0 / 0.70)

    def test_nothing_bound_gives_no_measured_value_and_falls_back_to_zero_load(self) -> None:
        """With no bound head there is no load to average, so the measured value must be NaN, not 0."""
        result = _bound_head_lifetime(_FakeSlice(np.zeros(8, np.int32), np.zeros(8)),
                                      _catch_slip_params())
        assert result["n_bound"] == 0
        assert math.isnan(float(result["mean_load_s"]))
        assert result["conservative_s"] == pytest.approx(result["zero_load_s"])

    def test_load_lengthens_a_catch_bond_and_the_conservative_value_follows(self) -> None:
        """A CATCH bond under load detaches more slowly, so the requirement can only get longer."""
        bound = np.ones(64, np.int32)
        loaded = _bound_head_lifetime(_FakeSlice(bound, np.full(64, 0.5)), _catch_slip_params())
        unloaded = _bound_head_lifetime(_FakeSlice(bound, np.zeros(64)), _catch_slip_params())
        assert loaded["mean_load_s"] > unloaded["mean_load_s"]
        assert loaded["conservative_s"] == pytest.approx(loaded["mean_load_s"])

    def test_load_shortens_a_slip_bond_so_the_zero_load_value_stays_conservative(self) -> None:
        """A pure Bell slip detaches FASTER under load, so the pre-declarable zero-load value governs."""
        bound = np.ones(64, np.int32)
        result = _bound_head_lifetime(_FakeSlice(bound, np.full(64, 2.0)), _slip_params())
        assert result["mean_load_s"] < result["zero_load_s"]
        assert result["conservative_s"] == pytest.approx(result["zero_load_s"]), (
            "taking the shorter of the two would let a loaded run declare settlement sooner, which is "
            "gate-loosening by arithmetic")

    def test_only_bound_heads_enter_the_average(self) -> None:
        """An unbound head has no bond to detach; including its load would bias the lifetime."""
        bound = np.array([1, 0, 1, 0], np.int32)
        loads = np.array([0.5, 99.0, 0.5, 99.0])
        mixed = _bound_head_lifetime(_FakeSlice(bound, loads), _catch_slip_params())
        only_bound = _bound_head_lifetime(_FakeSlice(np.ones(2, np.int32), np.array([0.5, 0.5])),
                                         _catch_slip_params())
        assert mixed["n_bound"] == 2
        assert mixed["mean_load_s"] == pytest.approx(only_bound["mean_load_s"])


class TestAssessSeries:
    """Every recorded observable is judged, and the SAMPLE spacing is what the estimator is given."""

    def test_all_observables_get_a_verdict(self) -> None:
        gamma = 4.0 + np.random.default_rng(0).normal(0.0, 0.05, size=4_000)
        reports = _assess_series(_trajectory(gamma), 0.01, driving_tau_s=1.0 / 0.7,
                                 min_windows=5.0, drift_sigma=1.0)
        assert set(reports) == set(STATIONARITY_OBSERVABLES)
        assert reports["gamma_total_pn_per_um"]["verdict"] == StationarityVerdict.STATIONARY.value
        assert reports["gamma_total_pn_per_um"]["mean"] == pytest.approx(4.0, abs=0.01)

    def test_the_2026_07_28_sweep_shape_is_refused(self) -> None:
        """60 samples of a still-climbing gamma at dt 0.01 is 0.6 s against a 1.43 s lifetime.

        Every point of that sweep must come back without a mean.  This is the regression test for the
        run that produced a non-monotonic dose response and was reported before it was judged.
        """
        climbing = np.linspace(1.5, 4.32, 60)
        reports = _assess_series(_trajectory(climbing), 0.01, driving_tau_s=1.0 / 0.7,
                                 min_windows=5.0, drift_sigma=1.0)
        for name in STATIONARITY_OBSERVABLES:
            assert reports[name]["mean"] is None, f"{name} returned a mean for a 0.6 s transient"
        assert reports["gamma_total_pn_per_um"]["verdict"] == StationarityVerdict.TOO_SHORT.value

    def test_a_telemetry_stride_changes_the_sample_spacing_not_the_step(self) -> None:
        """Feeding ``dt_phys`` instead of ``stride*dt_phys`` understates tau by exactly the stride.

        The stride exists because the per-step gamma readback dominates this driver's wall-clock at
        native; it must not buy that speed by quietly shrinking the measured correlation time.
        """
        rng = np.random.default_rng(1)
        phi = 0.9
        series = np.empty(20_000)
        series[0] = rng.normal()
        for k in range(1, series.size):
            series[k] = phi * series[k - 1] + rng.normal()
        traj = _trajectory(4.0 + 0.05 * series)
        correct = _assess_series(traj, 0.05, driving_tau_s=1.0 / 0.7, min_windows=5.0,
                                 drift_sigma=1.0)["gamma_total_pn_per_um"]
        wrong = _assess_series(traj, 0.01, driving_tau_s=1.0 / 0.7, min_windows=5.0,
                               drift_sigma=1.0)["gamma_total_pn_per_um"]
        assert correct["tau_int_s"] == pytest.approx(5.0 * wrong["tau_int_s"], rel=0.05)

    def test_too_few_samples_is_reported_not_raised(self) -> None:
        """A run stopped before four samples must record a refusal, not crash the artifact write."""
        reports = _assess_series(_trajectory(np.array([1.0, 2.0])), 0.01, driving_tau_s=1.0 / 0.7,
                                 min_windows=5.0, drift_sigma=1.0)
        assert reports["gamma_total_pn_per_um"]["verdict"] == StationarityVerdict.TOO_SHORT.value
        assert reports["gamma_total_pn_per_um"]["mean"] is None

    def test_reports_are_json_able(self) -> None:
        """These go straight into the run-record measurements block."""
        import json

        reports = _assess_series(_trajectory(np.linspace(0.0, 1.0, 500)), 0.01,
                                 driving_tau_s=1.0 / 0.7, min_windows=5.0, drift_sigma=1.0)
        json.dumps(reports, allow_nan=True, default=str)


def test_the_arithmetic_is_shared_not_driver_local() -> None:
    """The driver wrapper must delegate — a second lane has to inherit the fix, not re-derive it.

    This is the regression test for the structural defect, not for a number: ``_bound_head_lifetime``
    exists only to read the bound heads' loads off the device, and every value it reports must come from
    the shared closed form.
    """
    bound = np.ones(32, np.int32)
    loads = np.linspace(0.0, 1.5, 32)
    wrapper = _bound_head_lifetime(_FakeSlice(bound, loads), _catch_slip_params())
    shared = bound_state_lifetime(
        loads, catch_slip=True,
        k_catch0=CATCH_SLIP_PROVISIONAL["k_catch0"], x_catch=CATCH_SLIP_PROVISIONAL["x_catch"],
        k_slip0=CATCH_SLIP_PROVISIONAL["k_slip0"], x_slip=CATCH_SLIP_PROVISIONAL["x_slip"],
        kT=_catch_slip_params().kT)
    assert wrapper == shared


def test_the_shared_form_refuses_to_default_a_missing_rate() -> None:
    """A silently defaulted off-rate is the failure the shared function exists to prevent."""
    with pytest.raises(ValueError, match="refusing to default"):
        bound_state_lifetime(np.zeros(4), catch_slip=True, k_catch0=0.35, x_catch=1e-3, k_slip0=None,
                             x_slip=0.6e-3)
    with pytest.raises(ValueError, match="refusing to default"):
        bound_state_lifetime(np.zeros(4), catch_slip=False, k_off0=0.35)
