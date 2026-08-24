"""γ on the resting path — the properties that make it usable for `STATE.md` (e) 1, not its magnitude.

Nothing here asserts a γ value against a band. `STATE.md` (c) 3 and (c) 17 stand, and a test that pinned
a number would be the first place a retracted magnitude came back. What is pinned is the machinery (e) 1
needs to become decidable: the extraction is the ported one term for term, the instrument behaves like a
force per unit cut length, the window W is reported even when the answer is refused, and the error bar is
the SEED SCATTER — the one substitution `PI_DECISION_e1_STATIONARITY_2026-08-16.md` §3 D-2 calls the
single most likely way to write the amendment and still be measuring nothing.
"""

from __future__ import annotations

import ast
import pathlib
import subprocess
import sys

import numpy as np
import pytest
import warp as wp

from aleph.laws.gamma_estimator import method_of_planes_gamma
from aleph.observe.stationarity import (
    StationarityVerdict,
    assess_stationarity,
    integrated_autocorrelation_time,
)
from aleph.world import observe_gamma as mod
from aleph.world.observe_gamma import (
    NO_MAGNITUDE_CLAIMED,
    STEP_ACCEPTANCE_NOTE,
    actin_axial_tension,
    gamma_from_resting_readback,
    observe_gamma_seed_scatter,
    plane_force_sums,
)

DT = 0.01

CUDA = wp.get_cuda_device_count() > 0
cuda_only = pytest.mark.skipif(not CUDA, reason="requires a CUDA device")


def _ring(n: int = 48, radius: float = 5.0) -> tuple[np.ndarray, np.ndarray]:
    """A closed cortical ring of `n` nodes under a uniform inward pull — a cell-shaped test case."""
    theta = np.linspace(0.0, 2.0 * np.pi, n + 1)[:-1]
    pos = radius * np.stack([np.cos(theta), np.sin(theta), np.zeros_like(theta)], axis=1)
    inward = -pos / np.linalg.norm(pos, axis=1, keepdims=True)
    return pos, inward


def _settled(seed: int, level: float, *, n: int = 6000, rho: float = 0.97,
             sigma: float = 0.02, drift: float = 0.0) -> np.ndarray:
    """A correlated series settled tightly about `level`, optionally still going somewhere."""
    rng = np.random.default_rng(seed)
    out = np.zeros(n)
    for i in range(1, n):
        out[i] = rho * out[i - 1] + rng.standard_normal()
    return level + sigma * out + drift * np.linspace(0.0, 1.0, n)


def test_self_check_passes() -> None:
    mod._demo()


def test_actin_extraction_matches_the_ported_loop() -> None:
    """The vectorised running sum IS `assembled_cortical_stress.actin_axial_tension`, term for term.

    `world/` may not import `components/`, so the arithmetic had to be re-expressed rather than called.
    That is exactly the move that silently drifts, which is why the loop is written out here as the
    reference and compared, on a ragged case that includes a one-node filament and an empty one.
    """
    rng = np.random.default_rng(20260821)
    off = np.array([0, 7, 7, 8, 19, 40])          # 7-node, EMPTY, 1-node, 11-node, 21-node
    pos, f = rng.normal(size=(40, 3)), rng.normal(size=(40, 3))

    ref_a, ref_t = [], []
    for i in range(off.size - 1):
        s, e = int(off[i]), int(off[i + 1])
        if e - s < 2:
            continue
        cum = np.zeros(3)
        for k in range(e - s - 1):
            cum = cum + f[s + k]
            u = pos[s + k + 1] - pos[s + k]
            ref_a.append(s + k)
            ref_t.append(-float(np.dot(cum, u / np.linalg.norm(u))))

    seg_a, seg_b, tension = actin_axial_tension(f, pos, off)
    assert seg_a.tolist() == ref_a
    assert (seg_b == seg_a + 1).all()
    # Relative, not bit-exact: the loop accumulates per filament and this accumulates once and
    # subtracts, so the two reassociate the same float64 sum differently.
    assert np.allclose(tension, ref_t, rtol=1e-11, atol=1e-11)


def test_gamma_is_a_force_per_unit_cut_length() -> None:
    """Same elements, twice the cut circumference, half the γ. Pins the instrument's dimension."""
    pos, inward = _ring()
    kw = dict(pos=pos, f_ext=inward, n_actin=pos.shape[0],
              fiber_offsets=[0, pos.shape[0]], force_mask=())
    one = gamma_from_resting_readback(R_um=5.0, **kw)
    two = gamma_from_resting_readback(R_um=10.0, **kw)
    assert one.gamma_total_pn_per_um > 0.0
    assert two.gamma_total_pn_per_um == pytest.approx(0.5 * one.gamma_total_pn_per_um, rel=1e-12)


def test_the_cut_centre_does_not_depend_on_where_actin_sits_in_the_arena() -> None:
    """γ must not change when the SAME cortex is claimed at a different arena offset.

    ⚠ The regression this pins cost a driver A/B. The default centre was `pos[:n_actin].mean(axis=0)`,
    which is the actin centroid only if actin starts at index 0. The arena claims a cortex at
    `[lo, lo+N)` and addresses it by ABSOLUTE offsets, so a correct caller passes `lo+N` as `n_actin` —
    and that slice then averaged every membrane and envelope node below `lo` into the "actin centroid".
    A centre offset does not perturb the sum; it FLIPS WHICH ELEMENTS CROSS a plane, which is why it
    showed up as 7.5e-06 at step 0 against 3.865e-12 for the same code on identical inputs, growing to
    15% as the cell moved.
    """
    rng = np.random.default_rng(31)
    n, lo = 60, 200
    theta = np.linspace(0.0, 2.0 * np.pi, n + 1)[:-1]
    cortex = 5.0 * np.stack([np.cos(theta), np.sin(theta), np.zeros(n)], axis=1)
    f_cortex = -cortex / np.linalg.norm(cortex, axis=1, keepdims=True)

    at_zero = gamma_from_resting_readback(
        R_um=5.0, pos=cortex, f_ext=f_cortex, n_actin=n, fiber_offsets=[0, n], force_mask=())

    # The same cortex, claimed at offset `lo`, with foreign nodes (membrane, envelope) below it — placed
    # deliberately off-centre, which is what made the contaminated centroid move.
    pos = np.vstack([rng.normal(size=(lo, 3)) * 3.0 + np.array([8.0, -4.0, 2.0]), cortex])
    f_ext = np.vstack([rng.normal(size=(lo, 3)), f_cortex])
    at_offset = gamma_from_resting_readback(
        R_um=5.0, pos=pos, f_ext=f_ext, n_actin=lo + n, fiber_offsets=[lo, lo + n], force_mask=())

    assert at_offset.gamma_total_pn_per_um == pytest.approx(
        at_zero.gamma_total_pn_per_um, rel=1e-12), (
        "the same cortex measured at a different arena offset must give the same gamma")
    assert at_offset.centre_um == pytest.approx(at_zero.centre_um, abs=1e-12)
    # And the Newton-closure diagnostic must not have swept the foreign nodes in either.
    assert at_offset.net_force_vector_pn == pytest.approx(at_zero.net_force_vector_pn, abs=1e-12)


def test_an_absent_family_reports_none_not_zero() -> None:
    """A group with no elements must not report 0.0 — absent and measured-zero are different findings.

    ⚠ The PHASE 4 cell the 2026-08-21 tau runs were measured on stands 11 populations and never calls
    `build_nmii`: there is no motor in it. A `gamma_source` of 0.0 there says "the motors contributed
    nothing" when the truth is "there are no motors", and a reader would need the population list by
    heart to tell them apart. `None` cannot be summed, averaged or plotted as a zero.
    """
    pos, inward = _ring()
    kw = dict(R_um=5.0, pos=pos, f_ext=inward, n_actin=pos.shape[0],
              fiber_offsets=[0, pos.shape[0]], force_mask=())
    absent = gamma_from_resting_readback(**kw)
    assert absent.gamma_source_pn_per_um is None and absent.n_source_elements == 0
    assert absent.gamma_network_pn_per_um is not None
    assert absent.as_dict()["gamma_source_pn_per_um"] is None

    # A family that IS present and genuinely sums to zero must report a number, not None — otherwise
    # the distinction runs the other way and a real measurement is hidden.
    n = pos.shape[0]
    rA, rB = pos, np.roll(pos, -1, axis=0)
    present = gamma_from_resting_readback(
        **kw, extra_families={"crossbridge": (rA, rB, np.zeros(n))})
    assert present.gamma_source_pn_per_um == 0.0, "present-and-zero is a measurement"
    assert present.n_source_elements == n


def test_the_force_mask_cannot_be_defaulted() -> None:
    """A measurement masked differently from the solve reads a different force field (2026-07-29)."""
    pos, inward = _ring()
    with pytest.raises(TypeError, match="force_mask"):
        gamma_from_resting_readback(  # type: ignore[call-arg]
            R_um=5.0, pos=pos, f_ext=inward, n_actin=pos.shape[0], fiber_offsets=[0, pos.shape[0]])


def test_no_magnitude_claimed_survives_serialisation() -> None:
    """The stamp has to be in the JSON, not only in the docstring — the record is what gets quoted."""
    pos, inward = _ring()
    gamma = gamma_from_resting_readback(R_um=5.0, pos=pos, f_ext=inward, n_actin=pos.shape[0],
                                        fiber_offsets=[0, pos.shape[0]], force_mask=("steric_wca",))
    scatter = observe_gamma_seed_scatter(
        {s: _settled(s, 100.0 + 0.5 * i) for i, s in enumerate((1, 2, 3))}, DT)
    for record in (gamma.as_dict(), scatter.as_dict()):
        assert record["magnitude_claim"] == NO_MAGNITUDE_CLAIMED
        # `world/step.py` stands on `UndefinedAcceptance`; a series that did not say so would read as
        # sampled from accepted physics. The note travels with the record, not with the docstring.
        assert record["step_acceptance"] == STEP_ACCEPTANCE_NOTE
    assert gamma.as_dict()["force_mask"] == ["steric_wca"]


def test_the_error_bar_is_the_seed_scatter_and_the_within_run_sem_is_only_shown() -> None:
    """The whole point. `STATE.md` (f): compare γ against the seed scatter, never the within-run sem.

    The three replicates each settle tightly about their own level, which is the measured regime — so a
    within-run sem is a small, confident error bar around a number that moves when the seed does.
    """
    seeds = (11, 12, 13)
    series = {s: _settled(s, 100.0 + 0.5 * i) for i, s in enumerate(seeds)}
    got = observe_gamma_seed_scatter(series, DT)

    assert got.status == "MEASURED", got.refusals
    means = np.array([r["mean"] for r in got.per_replicate])
    assert got.mean_pn_per_um == pytest.approx(float(means.mean()))
    assert got.seed_scatter_pn_per_um == pytest.approx(float(means.std(ddof=1)))
    assert got.sem_across_seeds_pn_per_um == pytest.approx(
        float(means.std(ddof=1)) / np.sqrt(len(seeds)))

    # The failure mode, measured rather than asserted: the standing rule's factor is ~67.
    assert got.scatter_inflation > 10.0
    assert got.sem_across_seeds_pn_per_um > got.mean_within_run_sem_pn_per_um
    assert got.resolvable_fraction_one_seed > got.resolvable_fraction_n_seeds > 0.0
    assert "SEED SCATTER" in got.contract["error_bar_basis"]


def test_the_window_is_measured_and_reported_even_when_the_report_refuses() -> None:
    """(e) 1 is held on "stationarity of WHAT, over WHAT window". W is answerable either way."""
    two = observe_gamma_seed_scatter({1: _settled(1, 100.0), 2: _settled(2, 100.5)}, DT)
    assert two.status == "REFUSED" and two.mean_pn_per_um is None
    assert two.n_replicates == 2 and "minimum of 3" in " ".join(two.refusals)
    assert two.window_s is not None and two.window_s > 0.0
    assert two.equilibration_s is not None
    assert two.max_tau_int_s is not None and two.min_n_effective is not None
    assert two.all_sokal_windows_closed is not None


def test_a_drifting_replicate_refuses_the_report_rather_than_being_dropped() -> None:
    """Averaging the replicates that settled and discarding the one that did not is selection."""
    series = {s: _settled(s, 100.0 + 0.5 * i) for i, s in enumerate((21, 22, 23))}
    series[24] = _settled(24, 100.0, drift=50.0)
    got = observe_gamma_seed_scatter(series, DT)
    assert got.status == "REFUSED" and got.mean_pn_per_um is None
    assert got.n_replicates == 4, "the refusing replicate stays in the census"
    assert any("seed 24" in reason for reason in got.refusals), got.refusals


def test_a_constant_series_is_refused_however_long_it_gets() -> None:
    """The positive control for the worst failure this module can have: a perfect steady state of nothing.

    Measured 2026-08-21 on the resting path with `relax = 0`: gamma is constant to float64 round-off, and
    `assess_stationarity` on its own returns STATIONARY with `mean = 0.18860334` and `sem = 0.0`. At 20
    samples the length bar refuses it; at 300 it does not. So the refusal cannot rest on length, and this
    test asserts it at a length where the underlying verdict is STATIONARY — the case that would ship.
    """
    exact = np.full(6000, 0.18860334)
    solo = assess_stationarity(exact, DT, observable="control")
    assert solo.verdict is StationarityVerdict.STATIONARY and solo.sem == 0.0, (
        "if this ever stops being true the control below has stopped controlling for anything")

    # Constant to round-off rather than exactly: `zero_variance` is False, so the flag alone is not enough.
    rng = np.random.default_rng(0)
    near = 0.18858987 + rng.choice([-1, 0, 1], 6000) * np.spacing(0.18858987)
    assert assess_stationarity(near, DT, observable="control").verdict is StationarityVerdict.STATIONARY

    got = observe_gamma_seed_scatter({11: near, 12: exact, 13: exact + 1e-17}, DT)
    assert got.status == "REFUSED", "a constant gamma must never yield a mean"
    assert got.mean_pn_per_um is None and got.seed_scatter_pn_per_um is None
    assert got.degenerate_seeds == (11, 12, 13)
    assert all("MORE SAMPLES WILL NOT HELP" in r for r in got.refusals)


def test_the_window_is_also_reported_as_a_sample_count() -> None:
    """Seconds are `samples * dt_s`. When the caller's time axis is unsourced, the count still is not."""
    series = {s: _settled(s, 100.0 + 0.5 * i, n=4000) for i, s in enumerate((31, 32, 33))}
    got = observe_gamma_seed_scatter(series, DT)
    assert got.window_samples is not None and 0 < got.window_samples <= 4000
    assert got.window_s == pytest.approx(got.window_samples * DT)
    assert "only as sourced as" in got.contract["time_axis"]


def test_a_length_bar_that_would_admit_an_unclosed_sokal_window_is_refused() -> None:
    """The guarantee `all_sokal_windows_closed` can never be False in a MEASURED report, made structural.

    It holds at the defaults by arithmetic — an unclosed window bounds `n_eff` under `(2/3)*sokal_c` = 3.34
    against a required 25 — but that is a coincidence between two independently declared thresholds. A
    caller who lowers the length bar reopens the hole silently, so the contract is refused instead. This
    is not a threshold of ours: `4/3` falls out of the bound at the 4-sample minimum.
    """
    series = {s: _settled(s, 100.0 + 0.5 * i) for i, s in enumerate((51, 52, 53))}
    with pytest.raises(ValueError, match="unclosed Sokal window can pass"):
        observe_gamma_seed_scatter(series, DT, min_tau_windows=5.0, sokal_c=5.0)
    # Exactly at the bound is legal; the default (50 vs 6.67) clears it with room.
    observe_gamma_seed_scatter(series, DT, min_tau_windows=20.0 / 3.0, sokal_c=5.0)
    assert observe_gamma_seed_scatter(series, DT).all_sokal_windows_closed is True


def test_a_too_short_refusal_says_how_much_longer_in_sourced_units() -> None:
    """A refusal that names the deficit is actionable; one that only names the verdict costs a GPU night.

    `required_window_samples` is `min_tau_windows * tau_samples`, and samples are a count — sourced even
    when the run's `dt` is not. `min_tau_reliability` rides beside it because a `tau` measured over few
    correlation times is biased LOW, which makes both `n_effective` and this target optimistic.
    """
    series = {s: _settled(s, 100.0 + 0.5 * i, n=400, rho=0.99) for i, s in enumerate((61, 62, 63))}
    got = observe_gamma_seed_scatter(series, DT)
    assert got.status == "REFUSED"
    assert got.required_window_samples > got.window_samples, "a refusal must not ask for less than it had"
    assert got.required_window_samples == pytest.approx(
        np.ceil(50.0 * max(r["tau_samples"] for r in got.per_replicate)))
    assert 0.0 < got.min_tau_reliability < 50.0, (
        "this fixture is deliberately short relative to its own correlation time, which is what makes "
        "the reliability flag mean something here")


def test_seeds_drifting_apart_are_distinguishable_from_seeds_drifting_together() -> None:
    """Both are `REFUSED`. Only the SIGNS separate "still equilibrating" from "going elsewhere".

    `observe/stationarity.py` keeps the slope signed on purpose; a summary that dropped it applied `abs`
    one level up. Seeds going opposite ways is an input to D-3 — initial-condition independence, "without
    this, stationary is compatible with stuck where it started" — and this test pins that the two cases
    are told apart, NOT that either one means anything about the cell.

    ⚠⚠ **PI RULING 18 (b) NARROWED THIS, AND THE NARROWING IS THE POINT OF THE BLOCK BELOW.** Until
    2026-08-22 both groups were `DRIFTING` on every replicate and this test asserted the two verdict
    LISTS were equal — "told apart only by the signs". Filtering contaminated candidates moves the
    window much later, so the retained series went from 6,000 samples to 520–2,027, and at 520 a real
    drift stops being resolvable:

        seed 21  drift +0.0186/s  6,000 samples -> DRIFTING     520 samples -> 1.33 sigma, TOO_SHORT
        seed 24  drift +0.0118/s                                520 samples -> DRIFTING

    **The drift did not go away; the window got too short to see it.** That is a POSITIVE diagnosis
    becoming an INCONCLUSIVE one — the same class of loss (b) was ruled to prevent, arriving through a
    different door, and it lands on D-3's input rather than D-2's. Recorded here rather than smoothed
    over: the assertion that the two verdict lists match is GONE because it is no longer true, and what
    replaced it pins the property the test was actually written for.
    """
    def with_drift(seed: int, slope: float) -> np.ndarray:
        return _settled(seed, 100.0) + np.linspace(0.0, slope, 6000)

    opposed = observe_gamma_seed_scatter(
        {21: with_drift(21, +3.0), 22: with_drift(22, -3.0), 23: with_drift(23, +3.0)}, DT)
    together = observe_gamma_seed_scatter(
        {24: with_drift(24, 3.0), 25: with_drift(25, 3.0), 26: with_drift(26, 3.0)}, DT)

    assert opposed.status == together.status == "REFUSED"
    assert all(r["verdict"] != "MEASURED" for r in opposed.per_replicate + together.per_replicate), (
        "neither group may be MEASURED — whatever the verdict is called, a drifting series is not a "
        "measurement, and that is what survives 18 (b) here")

    # ── the property this test exists for, unchanged ──────────────────────────────────────────────
    assert opposed.drift_directions_agree is False
    assert together.drift_directions_agree is True
    assert opposed.drift_per_s_by_seed[21] > 0 > opposed.drift_per_s_by_seed[22]
    assert set(opposed.drift_sigma_by_seed) == {21, 22, 23}

    # ── and the narrowing, pinned so it cannot quietly deepen ─────────────────────────────────────
    # The SIGNIFICANCE is what the verdict is thresholded on; the slope carries only direction and size.
    # Before 18 (b) every seed here cleared 3 sigma. Now one does not, and the field that explains why
    # is the window length — so both are asserted together, because either alone reads as noise.
    assert opposed.drift_sigma_by_seed[21] < 3.0 < opposed.drift_sigma_by_seed[22], (
        "seed 21's drift stopped being resolvable after the filter shortened its window; if this ever "
        "changes, the measured block in the docstring above is stale and must be re-measured")
    assert opposed.window_samples is not None and opposed.window_samples < 1000, (
        "the cause is window length, not the drift — pin it beside the consequence")
    assert opposed.drift_per_s_by_seed[21] > 0.0, "the drift is still THERE; only its resolution went"
    assert opposed.as_dict()["drift_sigma_by_seed"][21] == opposed.drift_sigma_by_seed[21]
    assert opposed.as_dict()["drift_per_s_by_seed"][22] < 0, "the sign must survive serialisation"


def test_tau_carries_its_own_error_bar() -> None:
    """"tau grew" decides nothing; "tau grew by more than its own error bar" does.

    Madras-Sokal: `sigma(tau)/tau ~ sqrt(2*(2M+1)/N)`. A longer series measures the same correlation time
    with a smaller relative error, which is the property that makes a growth-vs-error comparison mean
    something, so that is what is pinned here rather than a value.
    """
    short = observe_gamma_seed_scatter(
        {s: _settled(s, 100.0 + 0.5 * i, n=1000) for i, s in enumerate((71, 72, 73))}, DT)
    long = observe_gamma_seed_scatter(
        {s: _settled(s, 100.0 + 0.5 * i, n=8000) for i, s in enumerate((71, 72, 73))}, DT)
    assert 0.0 < long.max_tau_relative_error < short.max_tau_relative_error
    assert long.as_dict()["max_tau_relative_error"] == long.max_tau_relative_error

    # It is the Madras-Sokal form, not merely something that shrinks with N. `window_lag` is not in
    # `as_dict()`, so the reference is recomputed from the series rather than read back off the report —
    # a check reconstructed from the report's own output would only prove the report is self-consistent.
    per = []
    for i, seed in enumerate((71, 72, 73)):
        series = _settled(seed, 100.0 + 0.5 * i, n=1000)
        tau = integrated_autocorrelation_time(series, DT)
        n_analysed = short.per_replicate[i]["n_samples"]
        per.append(np.sqrt(2.0 * (2.0 * tau.window_lag + 1.0) / n_analysed))
    assert short.max_tau_relative_error == pytest.approx(max(per), rel=0.5), (
        "the reported error must track the Madras-Sokal form computed independently; the tolerance is "
        "loose because the report's tau is measured on the post-equilibration window and this reference "
        "is measured on the whole series, so the window lags differ"
    )


def test_every_field_in_the_read_order_has_a_count_twin() -> None:
    """A caller with no physical time axis passes `dt_s = 1.0`, and then `*_s` names lie about units.

    ⚠ This is the cut-centre defect in a different dress: a CORRECT VALUE OF THE WRONG QUANTITY, which
    passes everything that checks values. `max_tau_int_s: 188` from a `dt_s = 1.0` run is 188 SAMPLES and
    reads as 188 seconds to anyone opening the record later. Contract III is exactly that case — the
    driver's `dt_phys_s = 0.05` is documented as "this driver's clock, NOT a physiological quantity", and
    the mobility is unsourced too — so the whole read order must be legible without knowing `dt_s`.
    """
    series = {s: _settled(s, 100.0 + 0.5 * i, n=6000) for i, s in enumerate((81, 82, 83))}
    counts = observe_gamma_seed_scatter(series, 1.0)
    seconds = observe_gamma_seed_scatter(series, 0.05)

    for field in ("window_samples", "max_tau_samples", "equilibration_samples",
                  "required_window_samples", "min_n_effective", "min_tau_reliability",
                  "max_tau_relative_error"):
        assert getattr(counts, field) is not None, field
        assert getattr(counts, field) == pytest.approx(getattr(seconds, field)), (
            f"{field} must not depend on dt_s — it is a count or a ratio")

    # ⚠ And the RATE scales the OTHER way — same suffix, reciprocal trap.
    rate_c = counts.drift_per_s_by_seed[81]
    rate_s = seconds.drift_per_s_by_seed[81]
    assert rate_s == pytest.approx(rate_c / 0.05, rel=1e-9), (rate_c, rate_s)
    assert "OPPOSITE way" in counts.contract["time_axis"]

    # ...while the `*_s` fields do scale with dt_s, which is what makes their names dangerous.
    assert seconds.window_s == pytest.approx(counts.window_s * 0.05)
    assert seconds.max_tau_int_s == pytest.approx(counts.max_tau_int_s * 0.05)
    assert counts.window_s == pytest.approx(counts.window_samples)
    assert "COUNT despite" in counts.contract["time_axis"]


def test_a_window_that_opens_on_a_transient_invalidates_the_tau_fields_visibly() -> None:
    """The tau fields are meaningless when the retained window still begins on a rise — say so at the top.

    ⚠ Measured 2026-08-21 on contract III seed 1: the report carried `max_tau_samples` 1,280.2 and
    `min_tau_reliability` 12.01, both computed on a window the equilibration search had itself flagged at
    9.28 sigma, while tau on the clean part of the same series is ~102 — a factor 12.6. They were read at
    face value because nothing at this level said they had been invalidated. Per-replicate the flag
    existed; the summary dropped it, which is the same "computed but not said" shape as the drift signs.
    """
    # ⚠ THIS TEST'S PREMISE CHANGED ON 2026-08-22 AND THE ASSERTION MOVED WITH IT, rather than being
    # relaxed. It used to rely on `limit = x.size // 2`: a transient longer than half the series had no
    # reachable cut, so the search RETURNED a contaminated window and had to say so. PI ruling 18 (b)
    # both widened the bound and made the search SKIP contaminated candidates, so for a series that
    # settles at all there is now a clean cut and `windows_open_on_transient` is empty. The invariant
    # under test is unchanged — **contamination must be visible at the top level** — but the field that
    # carries it is now `windows_skipped_on_transient`. Both cases are pinned below, because a fix that
    # moves a warning is only a fix if the new place is checked as hard as the old one was.
    n, knee = 6000, 3800
    def contaminated(seed: int) -> np.ndarray:
        t = np.arange(n) / knee
        ramp = 400.0 * np.minimum(t, 1.0) ** 0.5
        return _settled(seed, 100.0) + ramp
    rising = {s: contaminated(s) for s in (91, 92, 93)}
    got = observe_gamma_seed_scatter(rising, DT)
    assert got.status == "REFUSED"
    assert got.windows_skipped_on_transient, (
        "a search that walked around a rise must say so at the top level; after 18 (b) this is where "
        "that fact lives and windows_open_on_transient no longer carries it")
    assert {seed for seed, _ in got.windows_skipped_on_transient} <= {91, 92, 93}
    assert all(n_skipped > 0 for _, n_skipped in got.windows_skipped_on_transient)
    assert got.as_dict()["windows_skipped_on_transient"] == [
        list(t) for t in got.windows_skipped_on_transient]

    # ── and the ORIGINAL assertion, on the case that still produces it ────────────────────────────
    # A series that never settles has no clean candidate at all, so the search falls back and the flag
    # means what it always meant. Without this the old protection would have been deleted rather than
    # moved, which is the failure mode this whole test exists to catch.
    # A knee this late leaves no candidate whose window is clean: `MIN_RETAINED_SAMPLES` caps the cut
    # at 5,488, and every tail from there back straddles the rise. Measured, not assumed —
    # `fell_back_to_contaminated` is True at knee 5,600 and False at 3,800.
    def barely_settles(seed: int) -> np.ndarray:
        t = np.arange(n) / 5600
        return _settled(seed, 100.0) + 400.0 * np.minimum(t, 1.0) ** 0.5
    always = observe_gamma_seed_scatter({s: barely_settles(s) for s in (97, 98, 99)}, DT)
    assert always.status == "REFUSED"
    assert always.windows_open_on_transient, (
        "when NO candidate is clean the search must fall back AND still raise the flag")
    assert any("OPENS ON A TRANSIENT" in r for r in always.refusals), always.refusals

    # ⚠ A settled series is NOT skip-free, and finding that out is why this block reads as it does.
    # Measured: calm skips (12, 9, 3) against the contaminated series' 28 each. Noise trips the
    # detector on a few early candidates whatever the shape, so **non-empty is not the signal — the
    # COUNT is**, and an assertion of `== ()` would have been a wrong bar that happened to be strict.
    calm = observe_gamma_seed_scatter(
        {s: _settled(s, 100.0 + 0.5 * i) for i, s in enumerate((94, 95, 96))}, DT)
    assert calm.windows_open_on_transient == (), "a settled series must never raise the flag itself"
    assert max((k for _, k in calm.windows_skipped_on_transient), default=0) < min(
        k for _, k in got.windows_skipped_on_transient), (
        "the skip COUNT must separate a long dirty region from ordinary noise; if these overlap the "
        "field cannot do the job it was added for")


def test_the_invalidation_flags_come_first_in_the_record() -> None:
    """A reader that truncates long fields must still see what says "do not read the rest".

    ⚠ On 2026-08-21 `opens_on_transient: True` and `transient_ratio: 9.28` reached a terminal as the
    string `per_replicate (4,712 chars)`, and the tau fields beside them were read at face value. Pinning
    individual keys to the front had already been tried there and failed by omission — two were pinned,
    the third was not thought of. So the flags are top-level keys AND first, and this pins the order.
    """
    series = {s: _settled(s, 100.0 + 0.5 * i) for i, s in enumerate((101, 102, 103))}
    keys = list(observe_gamma_seed_scatter(series, DT).as_dict())
    # `windows_skipped_on_transient` joined the front block 2026-08-22. It is not a new kind of field —
    # it is where `windows_open_on_transient`'s signal WENT under PI ruling 18 (b), so leaving it out of
    # the block would reproduce the exact defect this test was written for: the flag that invalidates
    # the tau fields reaching a reader as `per_replicate (4,712 chars)`.
    assert keys[:7] == ["status", "refusals", "windows_open_on_transient",
                        "windows_skipped_on_transient", "degenerate_seeds",
                        "magnitude_claim", "step_acceptance"], keys[:7]
    assert len(keys) == len(set(keys)), "a duplicated key would make the record reader-dependent"
    # Every tau-derived field must come AFTER the flag that can invalidate it.
    flag = keys.index("windows_open_on_transient")
    for field in ("max_tau_samples", "min_tau_reliability", "max_tau_relative_error",
                  "required_window_samples", "min_n_effective"):
        assert keys.index(field) > flag, field


def test_a_scatter_needs_at_least_two_numbers_to_exist() -> None:
    with pytest.raises(ValueError, match="at least 2"):
        observe_gamma_seed_scatter({1: _settled(1, 100.0)}, DT, min_replicates=1)


def _elements(n: int, seed: int = 5) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """`n` load-bearing elements on a shell, with one deliberately degenerate zero-length member."""
    rng = np.random.default_rng(seed)
    u = rng.normal(size=(n, 3))
    rA = 5.0 * u / np.linalg.norm(u, axis=1, keepdims=True)
    rB = rA + rng.normal(size=(n, 3)) * 0.4
    if n > 3:
        rB[3] = rA[3]
    tension = rng.normal(size=n)
    centre = 0.5 * (rA.mean(axis=0) + rB.mean(axis=0))
    return rA, rB, tension, centre


@pytest.mark.parametrize("n", [1, 5, 1000, 50_000])
def test_the_plane_sums_reproduce_the_shared_estimator(n: int) -> None:
    """`plane_force_sums` IS `laws.gamma_estimator.method_of_planes_gamma`, collapsed one step later.

    The module claims a resting γ and a GATE-B γ differ by the cell state and not by the instrument.
    Re-expressing the estimator to get the per-plane vector is exactly the move that could quietly break
    that claim, so it is asserted here against the shared function rather than argued in a docstring.
    Measured agreement is 0 to 1 ULP; the tolerance is set for the summation order, not for slack.
    """
    rA, rB, tension, centre = _elements(n)
    reference = method_of_planes_gamma(rA, rB, tension, 5.0, n_planes=64, centre=centre)
    collapsed = float(np.mean(np.abs(
        plane_force_sums(rA, rB, tension, R_um=5.0, centre=centre, n_planes=64))))
    assert collapsed == pytest.approx(reference, rel=1e-13, abs=1e-300)


def test_group_gamma_is_additive_over_families_but_the_collapsed_gamma_is_not() -> None:
    """The property the one-pass-per-family design rests on, and the trap next to it.

    The PLANE SUM is additive over families; ``|·|`` and the orientation mean are not. Adding per-family
    γ VALUES would be a different — and wrong — number, so the second assertion pins that they differ
    and that the split version is the one matching a single estimator call on the concatenation.
    """
    rA, rB, tension, centre = _elements(50_000)
    kw = dict(R_um=5.0, centre=centre, n_planes=64)
    left = plane_force_sums(rA[:20_000], rB[:20_000], tension[:20_000], **kw)
    right = plane_force_sums(rA[20_000:], rB[20_000:], tension[20_000:], **kw)

    whole = method_of_planes_gamma(rA, rB, tension, 5.0, n_planes=64, centre=centre)
    assert float(np.mean(np.abs(left + right))) == pytest.approx(whole, rel=1e-13)

    summed_gammas = float(np.mean(np.abs(left))) + float(np.mean(np.abs(right)))
    assert summed_gammas != pytest.approx(whole, rel=1e-3), (
        "adding per-family gamma VALUES must not accidentally agree — if it does this test has stopped "
        "guarding the distinction it exists for"
    )


def test_the_chunk_size_is_performance_only() -> None:
    """A chunked reduction that changed the answer would be a silent, load-dependent measurement."""
    rA, rB, tension, centre = _elements(20_000)
    kw = dict(R_um=5.0, centre=centre, n_planes=64)
    one_block = plane_force_sums(rA, rB, tension, chunk=10**9, **kw)
    for chunk in (1, 7, 4096):
        assert np.allclose(plane_force_sums(rA, rB, tension, chunk=chunk, **kw),
                           one_block, rtol=1e-13, atol=1e-300)


def test_an_empty_family_contributes_the_additive_identity() -> None:
    """Empty is not a special case — it is zeros, which is what "adds nothing" has to mean here."""
    empty = plane_force_sums(np.zeros((0, 3)), np.zeros((0, 3)), np.zeros(0),
                             R_um=5.0, centre=np.zeros(3), n_planes=64)
    assert empty.shape == (64,) and not empty.any()


# ---------------------------------------------------------------------------------------------
# The device path. Everything below the type-check needs a GPU and is skipped without one.
# ---------------------------------------------------------------------------------------------

_KERNELS = (
    "_actin_axial_tension_kernel",
    "_crosslink_tension_kernel",
    "_plane_partial_sums_kernel",
    "_plane_reduce_kernel",
    "_endpoint_sum_partial_kernel",
)


def test_the_kernels_type_check_without_a_device() -> None:
    """Warp's own forward type inference over each kernel body. No launch, no CPU binary, no number.

    This is the ONLY verification available on the dev Mac and it is worth having precisely because the
    numeric gate below cannot run here: a type or shape regression in a kernel nobody can execute locally
    would otherwise surface hours later on the GPU host, inside someone's allocation.
    """
    for name in _KERNELS:
        kernel = getattr(mod, name)
        kernel.adj.build(None)
        assert kernel.adj.blocks and kernel.adj.blocks[0].body_forward, name


def test_the_device_path_refuses_a_non_cuda_device() -> None:
    """The charter has no CPU simulation path, so this must raise rather than quietly become one."""
    if CUDA:
        pytest.skip("this machine has CUDA; the refusal is only reachable without one")
    with pytest.raises(RuntimeError, match="non-CUDA device"):
        mod.gamma_planes_device(
            pos_d=None, f_ext_d=None, fiber_offsets_d=None, seg_offsets_d=None,
            n_segments=0, n_actin=0, R_um=5.0, centre=[0.0, 0.0, 0.0])


def test_segment_offsets_agrees_with_the_host_extraction() -> None:
    """The kernel's output ranges must match, fiber for fiber, what the host extraction produces."""
    rng = np.random.default_rng(11)
    off = np.array([0, 7, 7, 8, 19, 40])
    pos, f = rng.normal(size=(40, 3)), rng.normal(size=(40, 3))
    seg_off = mod.segment_offsets(off)
    assert seg_off[0] == 0 and len(seg_off) == len(off)
    assert seg_off[-1] == actin_axial_tension(f, pos, off)[0].size
    # A fiber with fewer than two nodes owns an EMPTY range rather than being dropped, so the thread
    # index stays the fiber index and no compaction map is needed.
    assert seg_off[2] == seg_off[1]


def test_nothing_is_defined_after_the_main_block() -> None:
    """A definition below `if __name__ == "__main__"` is invisible to `python -m`, and only to that.

    ⚠ This is a STATIC check because the dynamic one does not work. The device kernels,
    `segment_offsets` and `gamma_planes_device` were appended after the `__main__` block, so `python -m
    aleph.world.observe_gamma` called `_demo()` before any of them was bound. It died with a `NameError`
    in the production runtime — and it exits 0 here, because without CUDA `_demo` never enters the branch
    that touches those names. Reproduced both ways before writing this: the broken ordering still exits 0
    on the dev Mac.

    So neither the import-based test (module body finishes first, every name exists) nor a subprocess run
    on a machine without a card can see it. What can is the file's own shape, and that is checkable
    anywhere.
    """
    tree = ast.parse(pathlib.Path(mod.__file__).read_text())
    main = [n for n in tree.body if isinstance(n, ast.If) and "__main__" in ast.dump(n.test)]
    assert len(main) == 1, "one __main__ block expected"
    late = [
        getattr(n, "name", type(n).__name__)
        for n in tree.body
        if n.lineno > main[0].lineno and not isinstance(n, ast.If)
    ]
    assert not late, (
        f"defined after the __main__ block and therefore unbound when `python -m` runs it: {late}"
    )


def test_the_documented_invocation_actually_runs() -> None:
    """`python -m aleph.world.observe_gamma` as a SUBPROCESS, because that is how the module is run.

    It does NOT catch the ordering fault above on a machine without CUDA — see that test. What it does
    catch is anything that breaks the invocation on a path which executes here, and it asserts the module
    always says which half of its gates ran: silence about the device half is the failure mode the whole
    `skipped`-reads-as-*not applicable* finding was about.
    """
    done = subprocess.run(
        [sys.executable, "-m", "aleph.world.observe_gamma"],
        capture_output=True, text=True, timeout=600,
        cwd=str(pathlib.Path(__file__).resolve().parents[3]),
    )
    assert done.returncode == 0, f"stdout:\n{done.stdout}\nstderr:\n{done.stderr}"
    assert "self-check OK" in done.stdout, done.stdout
    assert ("DEVICE gates PASS" in done.stdout) or ("DEVICE GATES NOT RUN" in done.stdout), done.stdout


@cuda_only
def test_device_self_check_passes() -> None:
    """The device gates live in the MODULE, not here, and this delegates rather than copying them.

    `aleph/world/observe_gamma.py::_device_demo` holds the equivalence, chain-closure and bitwise
    reproducibility checks as plain asserts, because the production runtime (`ffn_sim`) has no pytest and
    no pytest process has ever run on the GPU host — so a gate that exists only as a CUDA-marked test is
    a gate that has never executed where the numbers are made. Keeping the arithmetic here too would put
    the authoritative copy in the place that cannot run it.
    """
    mod._device_demo()


def test_the_device_gate_refuses_to_run_without_cuda() -> None:
    """`_device_demo` must raise rather than quietly measure nothing when there is no card."""
    if CUDA:
        pytest.skip("this machine has CUDA; the refusal is only reachable without one")
    with pytest.raises(RuntimeError, match="no CPU path"):
        mod._device_demo()


def test_a_fully_transient_series_is_refused_but_NOT_by_the_transient_fields() -> None:
    """The contamination fields are not a sufficient read on their own, and here is the case that shows it.

    ⚠ Written after getting this wrong. PI ruling 18 (b) turned `window_opens_on_a_transient` from a
    LABEL on the winning window into a FILTER inside the selection loop, and that function has a
    declining branch — `if tenth / max(tau, _TAU_SAMPLES_FLOOR) < _MIN_TRANSIENT_BLOCK_TAUS:
    return (False, 0.0)`. On a series still rising at its OWN END the trailing tenth is itself a rise,
    its tau is enormous, and the detector declines for nearly every candidate. Measured 2026-08-22 on
    one shape at three knees:

        knee 3,800 -> start 3,836, 28 skipped, MEASURED          (filter working)
        knee 5,600 -> start 5,480, 41 skipped, REFUSED, flag up  (fallback working)
        knee 5,900 -> start     0, 14 skipped, REFUSED, flag DOWN

    The third keeps the whole transient and both contamination fields read almost clean. ⚠ **It is
    still REFUSED** — the drift test independently calls it at 13.47 sigma, which is exactly what the
    declining branch's own docstring promised would happen ("such a series is TOO_SHORT, and the length
    gate will say so"). So this is defence-in-depth working, NOT a false pass, and it is recorded as a
    passing test rather than an xfail ratchet for that reason.

    What it pins is the thing a reader can get wrong: **an empty `windows_open_on_transient` plus a
    small `windows_skipped_on_transient` does not mean settled.** Read the verdict, not the flags.
    """
    n = 6000
    def all_transient(seed: int) -> np.ndarray:
        t = np.arange(n) / 5900
        return _settled(seed, 100.0) + 400.0 * np.minimum(t, 1.0) ** 0.5
    got = observe_gamma_seed_scatter({s: all_transient(s) for s in (81, 82, 83)}, DT)

    assert got.status == "REFUSED", "a series that rises from end to end must never be MEASURED"
    assert got.equilibration_samples == 0, (
        "the search kept everything — that is the condition this test is about, and if it ever stops "
        "being true the case has changed and the docstring above is stale")
    assert got.windows_open_on_transient == (), "the flag is DOWN here; that is the point"
    assert any("DRIFTING" in r for r in got.refusals), (
        "the refusal must come from the drift test, because the transient fields did not carry it")
