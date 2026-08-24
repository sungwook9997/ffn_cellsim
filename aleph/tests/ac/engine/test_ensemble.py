"""Structural gates for the ensemble validation harness + event-log schema (spec §7).

Pure host, no ``warp``, no CUDA.  Covers: EventRecord/EventLog validation + round-trip; reproducible +
independent per-realization seeds; the caller-owned ``simulate`` aggregation into an EnsembleSummary
(reproducible when the spy uses the passed rng); ObservableDistribution band/quantile correctness;
mix_ensembles proportion discipline; and the runner's tag/index guard.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.engine import (
    EnsembleRunner,
    EnsembleSummary,
    EventLog,
    EventRecord,
    ObservableDistribution,
    RealizationResult,
    mcf7_reference_state,
    mix_ensembles,
)


# --------------------------------------------------------------------------------------------------
# EventRecord
# --------------------------------------------------------------------------------------------------
def test_event_record_position_xor_region_amplitude_xor_rate():
    rec = EventRecord(
        event_type="nucleation",
        start_time=0.5,
        lifetime=1.0,
        mechanistic_trigger="g_actin_binding",
        position=(1.0, 2.0, 3.0),
        amplitude=4.0,
    )
    assert rec.position == (1.0, 2.0, 3.0)
    assert rec.region is None
    assert rec.amplitude == 4.0 and rec.rate_modifier is None

    rec2 = EventRecord(
        event_type="myosin_pulse",
        start_time=0.0,
        lifetime=0.0,
        mechanistic_trigger="signaling",
        region="cortex_apical",
        rate_modifier=2.5,
    )
    assert rec2.region == "cortex_apical" and rec2.position is None
    assert rec2.rate_modifier == 2.5 and rec2.amplitude is None


def test_event_record_rejects_both_and_neither_location():
    with pytest.raises(ValueError):
        EventRecord("t", 0.0, 0.0, "trig", position=(0.0, 0.0, 0.0), region="r", amplitude=1.0)
    with pytest.raises(ValueError):
        EventRecord("t", 0.0, 0.0, "trig", amplitude=1.0)  # neither position nor region


def test_event_record_rejects_both_and_neither_magnitude():
    with pytest.raises(ValueError):
        EventRecord("t", 0.0, 0.0, "trig", region="r", amplitude=1.0, rate_modifier=2.0)
    with pytest.raises(ValueError):
        EventRecord("t", 0.0, 0.0, "trig", region="r")  # neither amplitude nor rate_modifier


def test_event_record_rejects_negative_and_nan_times():
    with pytest.raises(ValueError):
        EventRecord("t", -1.0, 0.0, "trig", region="r", amplitude=1.0)
    with pytest.raises(ValueError):
        EventRecord("t", 0.0, -0.1, "trig", region="r", amplitude=1.0)
    with pytest.raises(ValueError):
        EventRecord("t", float("nan"), 0.0, "trig", region="r", amplitude=1.0)
    with pytest.raises(ValueError):
        EventRecord("", 0.0, 0.0, "trig", region="r", amplitude=1.0)  # empty type


def test_event_record_to_dict_round_trips():
    rec = EventRecord(
        event_type="sever",
        start_time=1.5,
        lifetime=0.25,
        mechanistic_trigger="cofilin",
        position=(0.1, 0.2, 0.3),
        rate_modifier=0.7,
    )
    d = rec.to_dict()
    assert d == {
        "event_type": "sever",
        "start_time": 1.5,
        "lifetime": 0.25,
        "mechanistic_trigger": "cofilin",
        "position": [0.1, 0.2, 0.3],
        "region": None,
        "amplitude": None,
        "rate_modifier": 0.7,
    }
    import json

    assert json.loads(json.dumps(d)) == d


# --------------------------------------------------------------------------------------------------
# EventLog
# --------------------------------------------------------------------------------------------------
def test_event_log_append_and_to_dict():
    state = mcf7_reference_state()
    log = EventLog(cell_state_tag=state, base_seed=7, realization_index=2)
    assert log.cell_state_tag == state.axes
    rec = EventRecord("nucleation", 0.0, 1.0, "g_actin", region="cortex", amplitude=1.0)
    log.append(rec)
    assert log.records == (rec,)
    d = log.to_dict()
    assert d["base_seed"] == 7 and d["realization_index"] == 2
    assert d["cell_state_tag"] == list(state.axes)
    assert d["records"][0]["event_type"] == "nucleation"


def test_event_log_validates_tag_and_seed():
    with pytest.raises(ValueError):
        EventLog(cell_state_tag=("a", "b"), base_seed=0, realization_index=0)  # not 6 axes
    with pytest.raises(ValueError):
        EventLog(cell_state_tag=mcf7_reference_state(), base_seed=-1, realization_index=0)
    with pytest.raises(ValueError):
        EventLog(cell_state_tag=("a", "b", "", "d", "e", "f"), base_seed=0, realization_index=0)


# --------------------------------------------------------------------------------------------------
# EnsembleRunner.seeds
# --------------------------------------------------------------------------------------------------
def test_seeds_deterministic_distinct_and_base_sensitive():
    runner = EnsembleRunner(cell_state=mcf7_reference_state(), base_seed=1234, n_realizations=8)
    seeds_a = runner.seeds()
    seeds_b = runner.seeds()
    assert seeds_a == seeds_b  # reproducible in base_seed
    assert len(set(seeds_a)) == len(seeds_a)  # distinct / independent
    other = EnsembleRunner(cell_state=mcf7_reference_state(), base_seed=1235, n_realizations=8)
    assert other.seeds() != seeds_a  # different base_seed -> different list


# --------------------------------------------------------------------------------------------------
# EnsembleRunner.run
# --------------------------------------------------------------------------------------------------
def _spy_simulate(rng, cell_state, index):
    """A tiny deterministic spy that DRAWS from the passed rng (so reproducibility is actually tested)."""
    connectivity = float(rng.uniform(90.0, 100.0))
    gamma = float(rng.normal(1.0, 0.1))
    lengths = rng.normal(3.0, 0.5, size=5)
    return RealizationResult(
        cell_state_tag=cell_state,
        seed=int(rng.integers(0, 2**31)),
        realization_index=index,
        scalars={"connectivity_pct": connectivity, "gamma": gamma},
        distributions={"length_um": lengths},
    )


def test_run_aggregates_summary():
    runner = EnsembleRunner(cell_state=mcf7_reference_state(), base_seed=42, n_realizations=6)
    summary = runner.run(_spy_simulate)
    assert isinstance(summary, EnsembleSummary)
    assert summary.n_realizations == 6
    assert len(summary.results) == 6
    # one ObservableDistribution per scalar, each with count == N
    assert set(summary.scalar_distributions) == {"connectivity_pct", "gamma"}
    for dist in summary.scalar_distributions.values():
        assert dist.count == 6
    # pooled distribution pools every realization's samples (6 * 5)
    assert summary.pooled_distributions["length_um"].count == 30
    assert summary.seeds == runner.seeds()


def test_run_is_reproducible_in_base_seed():
    a = EnsembleRunner(cell_state=mcf7_reference_state(), base_seed=99, n_realizations=5).run(_spy_simulate)
    b = EnsembleRunner(cell_state=mcf7_reference_state(), base_seed=99, n_realizations=5).run(_spy_simulate)
    np.testing.assert_array_equal(
        a.scalar_distributions["gamma"].values, b.scalar_distributions["gamma"].values
    )
    np.testing.assert_array_equal(
        a.pooled_distributions["length_um"].values, b.pooled_distributions["length_um"].values
    )
    # a different base seed gives different numbers
    c = EnsembleRunner(cell_state=mcf7_reference_state(), base_seed=100, n_realizations=5).run(_spy_simulate)
    assert not np.array_equal(
        a.scalar_distributions["gamma"].values, c.scalar_distributions["gamma"].values
    )


def test_run_rejects_mismatched_tag():
    other_tag = ("HeLa", "epithelial", "S", "glass", "round", "hypotonic")

    def bad_tag(rng, cell_state, index):
        return RealizationResult(
            cell_state_tag=other_tag, seed=1, realization_index=index, scalars={"x": 1.0}
        )

    runner = EnsembleRunner(cell_state=mcf7_reference_state(), base_seed=1, n_realizations=2)
    with pytest.raises(ValueError):
        runner.run(bad_tag)


def test_run_rejects_mismatched_index():
    def bad_index(rng, cell_state, index):
        return RealizationResult(
            cell_state_tag=cell_state, seed=1, realization_index=index + 100, scalars={"x": 1.0}
        )

    runner = EnsembleRunner(cell_state=mcf7_reference_state(), base_seed=1, n_realizations=2)
    with pytest.raises(ValueError):
        runner.run(bad_index)


# --------------------------------------------------------------------------------------------------
# ObservableDistribution
# --------------------------------------------------------------------------------------------------
def test_observable_distribution_within_band_and_quantiles():
    values = np.arange(1.0, 11.0)  # 1..10
    dist = ObservableDistribution(name="x", values=values)
    assert dist.count == 10
    assert dist.mean == pytest.approx(5.5)
    # ddof=1 std of 1..10
    assert dist.std == pytest.approx(float(np.std(values, ddof=1)))
    # within [3, 7] inclusive -> {3,4,5,6,7} = 5/10
    assert dist.within_band(3.0, 7.0) == pytest.approx(0.5)
    lo, med, hi = dist.quantiles((0.0, 0.5, 1.0))
    assert lo == pytest.approx(1.0) and hi == pytest.approx(10.0)
    assert med == pytest.approx(np.quantile(values, 0.5))


def test_observable_distribution_edge_cases():
    single = ObservableDistribution(name="x", values=np.array([2.0]))
    assert single.std == 0.0  # ddof guard
    empty = ObservableDistribution(name="x", values=np.array([]))
    assert empty.within_band(0.0, 1.0) == 0.0
    with pytest.raises(ValueError):
        ObservableDistribution(name="x", values=np.array([1.0])).within_band(2.0, 1.0)  # lo > hi


# --------------------------------------------------------------------------------------------------
# mix_ensembles
# --------------------------------------------------------------------------------------------------
def _constant_summary(tag, value, n=10):
    """An EnsembleSummary whose 'obs' scalar is constant == value across n realizations."""
    results = tuple(
        RealizationResult(cell_state_tag=tag, seed=i, realization_index=i, scalars={"obs": value})
        for i in range(n)
    )
    scalar = {"obs": ObservableDistribution(name="obs", values=np.full(n, value))}
    return EnsembleSummary(
        cell_state_tag=tag,
        n_realizations=n,
        seeds=tuple(range(n)),
        scalar_distributions=scalar,
        pooled_distributions={},
        results=results,
    )


def test_mix_ensembles_proportions_must_sum_to_one():
    tag_a = ("MCF7", "hybrid_E_M", "G1", "collagen_spread", "polarized", "physiological_resting")
    tag_b = ("MCF7", "hybrid_E_M", "S", "collagen_spread", "polarized", "physiological_resting")
    a = _constant_summary(tag_a, 0.0)
    b = _constant_summary(tag_b, 1.0)
    with pytest.raises(ValueError):
        mix_ensembles([(a, 0.5), (b, 0.4)])  # sums to 0.9
    with pytest.raises(ValueError):
        mix_ensembles([(a, 0.0), (b, 1.0)])  # non-positive proportion


def test_mix_ensembles_50_50_composition():
    tag_a = ("MCF7", "hybrid_E_M", "G1", "collagen_spread", "polarized", "physiological_resting")
    tag_b = ("MCF7", "hybrid_E_M", "S", "collagen_spread", "polarized", "physiological_resting")
    a = _constant_summary(tag_a, 0.0, n=10)
    b = _constant_summary(tag_b, 1.0, n=10)
    mixture = mix_ensembles([(a, 0.5), (b, 0.5)])
    obs = mixture["obs"]
    # 50/50 over a 20-sample pool -> 10 zeros + 10 ones
    zeros = np.count_nonzero(obs.values == 0.0)
    ones = np.count_nonzero(obs.values == 1.0)
    assert zeros == 10 and ones == 10
    assert obs.mean == pytest.approx(0.5)


def test_mix_ensembles_only_shared_observables():
    tag_a = ("MCF7", "hybrid_E_M", "G1", "collagen_spread", "polarized", "physiological_resting")
    tag_b = ("MCF7", "hybrid_E_M", "S", "collagen_spread", "polarized", "physiological_resting")
    a = _constant_summary(tag_a, 0.0)
    b = _constant_summary(tag_b, 1.0)
    # give `a` an extra observable that b lacks -> must not appear in the mixture
    a_extra = EnsembleSummary(
        cell_state_tag=tag_a,
        n_realizations=10,
        seeds=tuple(range(10)),
        scalar_distributions={
            "obs": ObservableDistribution(name="obs", values=np.zeros(10)),
            "extra": ObservableDistribution(name="extra", values=np.ones(10)),
        },
        pooled_distributions={},
        results=a.results,
    )
    mixture = mix_ensembles([(a_extra, 0.5), (b, 0.5)])
    assert set(mixture) == {"obs"}
