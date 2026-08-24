"""Two `stationarity.py` exist. This pins exactly how far they may differ.

`aleph/observe/stationarity.py` (the observation stack: `fluctuation`, `kinematics`, `tether`) and
`aleph/engine/observe/stationarity.py` (the production GATE-B driver) both define
`integrated_autocorrelation_time` and `assess_stationarity`. Nothing in the tree says which is canonical.

Measured 2026-08-20 rather than read off the source, because the first reading was wrong:

* the STATISTIC is identical — `tau_int` agrees to a ratio of 1.00000 on white noise, AR(1) rho=0.9,
  AR(1) rho=0.99 and a drifting series;
* the VERDICT CONTRACTS were not, and `min_windows` 5.0 vs `min_tau_windows` 50.0 were NOT the same
  knob: the engine's counts multiples of an EXTERNALLY supplied driving time, the other's multiples of
  the MEASURED `tau_int`. At the 2.5 s NMII driving time the external bound asks for 12.5 s, while
  50 x tau ranges 0.08-29 s over AR(1) rho 0.5-0.999 — neither dominates.

**PARTIALLY RESOLVED 2026-08-20 by PI decision** ("일단 B로 진행"), and deliberately only partially:

* ADOPTED — the self-referential refusal. The engine module gained `min_tau_windows = 50.0`
  (`n_eff >= 25`), which it simply lacked. `aleph/observe/stationarity.py` names the external-bound-only
  design as a failure mode: *"A series can clear that bar and still be hopeless for its own correlation
  time."* Both bounds are now kept and both are recorded in the contract, exactly as that module
  recommends. Disagreement fell 30.7% -> 24.8%.
* HELD — `drift_sigma`. The adopted module is LOOSER there (3.0 vs 1.0), and the decision was taken on
  a comparison this session got wrong. Applying it measurably flipped cases to STATIONARY and moved the
  dominant disagreement from 74 to 89 — loosening made agreement WORSE. A gate threshold moved on a
  mistaken premise is what the charter forbids, so it was reverted and surfaced. Open to the PI in
  `docs/v2_audit/ENGINE_FORWARD_ACCEPTANCE_2026-08-20.md`.

**What this test is for.** Not to bless the duplication — to stop it growing while it waits. If the shared
statistic ever diverges, that is a silent change to what every recorded `tau_int` in the archive means, and
the first test below fails. If the differing defaults are changed, the second fails and whoever changed
them has to say so in a diff rather than in a run six weeks later.

Delete this file when the two modules are consolidated; that is the intended end state, not this test.
"""

from __future__ import annotations

import numpy as np
import pytest

import aleph.engine.observe.stationarity as engine_mod
import aleph.observe.stationarity as observe_mod

DT = 1e-3


def _ar1(rho: float, n: int, seed: int, drift: float = 0.0) -> np.ndarray:
    """A first-order autoregressive series with a known correlation time, optionally drifting."""
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal(n)
    out = np.zeros(n)
    for i in range(1, n):
        out[i] = rho * out[i - 1] + noise[i]
    return out + drift * np.linspace(0.0, 1.0, n)


@pytest.mark.parametrize(
    ("label", "series"),
    [
        ("white", np.random.default_rng(0).standard_normal(4000)),
        ("ar1_0.9", _ar1(0.9, 4000, seed=1)),
        ("ar1_0.99", _ar1(0.99, 4000, seed=2)),
        ("drifting", np.linspace(0.0, 1.0, 4000) + 0.05 * np.random.default_rng(3).standard_normal(4000)),
    ],
)
def test_the_shared_statistic_is_still_shared(label: str, series: np.ndarray) -> None:
    """`tau_int` must stay identical. A divergence silently redefines every recorded tau_int."""
    tau_engine = float(engine_mod.integrated_autocorrelation_time(series, DT)[0])
    tau_observe = float(observe_mod.integrated_autocorrelation_time(series, DT).tau_int_s)
    assert tau_observe > 0.0, f"{label}: observe module returned a non-positive tau_int"
    assert tau_engine == pytest.approx(tau_observe, rel=1e-12), (
        f"{label}: the two stationarity modules no longer compute the same integrated autocorrelation "
        f"time — engine {tau_engine!r} vs observe {tau_observe!r}. Every tau_int already recorded was "
        "measured under the assumption that they agree."
    )


def test_the_verdict_contracts_still_differ_exactly_as_recorded() -> None:
    """Pin the differing defaults, so a change to either has to be visible in a diff."""
    import inspect

    engine_defaults = {
        name: param.default
        for name, param in inspect.signature(engine_mod.assess_stationarity).parameters.items()
    }
    observe_defaults = {
        name: param.default
        for name, param in inspect.signature(observe_mod.assess_stationarity).parameters.items()
    }
    # ADOPTED: the engine now carries the self-referential bound as well as its external one.
    assert engine_defaults["min_windows"] == 5.0, "the external driving-time bound was dropped"
    assert engine_defaults["min_tau_windows"] == 50.0, (
        "the self-referential bound adopted 2026-08-20 is gone; n_eff >= 25 is no longer enforced"
    )
    assert observe_defaults["min_tau_windows"] == 50.0
    # HELD: still different, on purpose. If these ever match, the open decision was taken — say so.
    assert engine_defaults["drift_sigma"] == 1.0, (
        "drift_sigma moved. It was held at 1.0 because the adopted module is LOOSER here and the "
        "decision rested on a mistaken comparison; moving it loosens a production gate."
    )
    assert observe_defaults["drift_sigma"] == 3.0
    assert engine_defaults["driving_correlation_time_s"] is inspect.Parameter.empty, (
        "the engine module's driving correlation time is a REQUIRED argument; if it gained a default, "
        "the production driver's window contract changed"
    )


def test_the_verdicts_still_disagree_at_the_recorded_rate() -> None:
    """The consequence, made concrete: the two modules disagree on ~a third of series.

    Measured 2026-08-20 over 420 AR(1) series: 129 disagreements (30.7%) before the self-referential
    bound was adopted, **104 (24.8%) after**, with BOTH directions still present — which is why this
    cannot be settled by "take the conservative module". Note the ordering that made the decision: also
    raising `drift_sigma` to the adopted module's 3.0 gave 27.4%, i.e. WORSE than adopting the refusal
    alone. A smaller deterministic sweep is used here so the test stays fast.

    Asserted as a FACT about today's tree, not as desirable behaviour. When the modules are consolidated,
    delete this file — that is the intended end state, not this test.
    """
    disagreements = 0
    total = 0
    for rho in (0.5, 0.9, 0.99):
        for drift in (0.0, 1.0, 2.0, 4.0):
            for seed in range(5):
                series = _ar1(rho, 4000, seed=seed, drift=drift)
                total += 1
                tau = float(engine_mod.integrated_autocorrelation_time(series, DT)[0])
                engine_verdict = engine_mod.assess_stationarity(
                    series, DT, observable="x", driving_correlation_time_s=tau
                ).verdict
                observe_verdict = observe_mod.assess_stationarity(series, DT, observable="x").verdict
                if str(engine_verdict) != str(observe_verdict):
                    disagreements += 1
    assert disagreements > 0, (
        f"the two stationarity modules now agree on all {total} swept series. That is good news, not a "
        "failure — but the recorded finding is then stale: re-measure, update "
        "docs/v2_audit/ENGINE_FORWARD_ACCEPTANCE_2026-08-20.md, and delete this file."
    )


def test_drift_sigma_3_would_be_a_gate_that_cannot_fire() -> None:
    """Why `drift_sigma` was HELD at 1.0 rather than raised to the adopted module's 3.0.

    Measured 2026-08-20 when the PI asked whether 1.0 is too strict. It is not — and 3.0 is not merely
    looser, it is above the statistic's practical ceiling.

    `drift_over_std` is computed on the window the equilibration search KEEPS, and that search cuts the
    drifting head off. Over AR(1) rho=0.9 series carrying a real linear drift of 3x and 6x their own
    standard deviation, the retained window's ratio maxes out at **2.73** and never reaches 3.0. A
    threshold above what the statistic can attain is a gate that cannot fail — the same defect this
    project diagnosed in `balance_ok`, which accepts a step when two force channels cancel to within
    float64 roundoff.

    Under the null (genuinely stationary, burned-in AR(1)) the false-DRIFTING rate at 1.0 is 0-2.4%, and
    the power at 1.0 is 52% at a 2-sigma drift and 94% at 3-sigma. So 1.0 sits inside the reachable range
    and 3.0 sits outside it.

    ⚠ Scope: this does NOT say a 3.0 run calls everything stationary. The transient guard still fires on
    an obvious ramp. What dies is the DRIFT TEST itself, so the cases that flip are the intermediate ones
    — real drift with no sharp transient for the guard to catch.
    """
    rng_std_drift = 6.0
    ratios = []
    for seed in range(30):
        burn = 30000
        n = 4000
        rng = np.random.default_rng(seed)
        noise = rng.standard_normal(n + burn)
        raw = np.zeros(n + burn)
        for i in range(1, n + burn):
            raw[i] = 0.9 * raw[i - 1] + noise[i]
        x = raw[burn:]
        series = x + rng_std_drift * float(np.std(x)) * np.linspace(0.0, 1.0, n)
        report = engine_mod.assess_stationarity(
            series, DT, observable="x", driving_correlation_time_s=1e-9
        )
        if np.isfinite(report.drift_over_std):
            ratios.append(float(report.drift_over_std))

    assert ratios, "no series produced a finite drift ratio; the probe itself is broken"
    peak = max(ratios)
    assert peak < 3.0, (
        f"`drift_over_std` reached {peak:.2f} on a {rng_std_drift}x-std drift, so a threshold of 3.0 is "
        "now attainable and this finding is stale — re-measure the power table before quoting it."
    )
    assert peak > 1.0, (
        f"`drift_over_std` peaked at {peak:.2f}, below the 1.0 threshold in force. If a 6x-std drift "
        "cannot reach 1.0 either, then the drift test is inert at BOTH settings and the criterion — not "
        "the constant — is what needs fixing."
    )
