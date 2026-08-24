"""A swept ``k_xb`` has a validity window, and the repo already knew where it was.

``ac/motor/minifilament_topology.working_stroke_strain`` calls itself **"the k_xb MASTER-knob arbiter"** and
states in its own docstring that the physiological working stroke of 5-20 nm needs ``k_xb`` ~100-1000 pN/µm,
while ``k_xb = 1 pN/µm`` gives a strain larger than the whole minifilament so "the stall geometry collapses".

That function existed the whole time the 2026-07-29 sweep was designed, queued and run, and nothing consulted
it.  Measured afterwards, with ``f_stall`` = 0.5 pN/head:

    k_xb   1000 -> 0.5 nm   100 -> 5.0 nm   10 -> 50 nm   3 -> 167 nm   1 -> 500 nm

Only ``k_xb=100`` is inside the window.  gamma's exponent flips sign between 50 nm and 167 nm
(-0.08 -> +1.20), and at 500 nm the strain exceeds the head offset itself (200 nm).  So the sweep's
headline shape is a property of running the model outside its stated domain — which is a finding about the
experiment, not a prediction of the model.

The window is imported from that arbiter's docstring rather than chosen here, which is the only reason it
is not just another threshold a session invented after seeing the data.
"""

from __future__ import annotations

import pytest

from aleph.scripts.ac_bleb_dose_response import (
    WORKING_STROKE_NM,
    _outside_working_stroke,
    refuse_unless_admissible,
)


def _pt(k_xb: float, mean: float = 4.0, sem: float = 0.001) -> dict:
    return {
        "name": f"kxb_{k_xb:g}", "k_xb": k_xb, "build": "aaaaaaaa", "steps_run": 800,
        "closure_sha256": "", "stationarity": {
            "gamma_total_pn_per_um": {"verdict": "STATIONARY", "mean": mean, "sem": sem,
                                      "drift_over_std": 0.2},
            "bound_fraction": {"verdict": "STATIONARY", "mean": 0.98, "sem": 1e-4},
        },
    }


def test_the_one_physiological_point_of_the_measured_sweep_is_k_xb_100() -> None:
    """5.0 nm at k_xb=100 — the only dose in the swept set that simulates a myosin head."""
    inside = {p["k_xb"] for p in [_pt(k) for k in (1000, 100, 10, 3, 1)]} - {
        p["k_xb"] for p, _ in _outside_working_stroke([_pt(k) for k in (1000, 100, 10, 3, 1)])}
    assert inside == {100}


@pytest.mark.parametrize("k_xb,nm", [(1000, 0.5), (100, 5.0), (10, 50.0), (3, 500.0 / 3), (1, 500.0)])
def test_the_strain_each_dose_implies(k_xb: float, nm: float) -> None:
    """f_stall / k_xb, in nm — the arithmetic the refusal rests on, pinned so it cannot drift silently."""
    out = dict((p["k_xb"], strain) for p, strain in _outside_working_stroke([_pt(k_xb)]))
    if WORKING_STROKE_NM[0] <= nm <= WORKING_STROKE_NM[1]:
        assert out == {}
    else:
        assert out[k_xb] == pytest.approx(nm, rel=1e-6)


def test_too_stiff_and_too_soft_are_both_refused_and_distinguished() -> None:
    """Sub-nm and half-micron strokes are opposite failures; a reader must not have to work out which."""
    reasons = refuse_unless_admissible([_pt(1000), _pt(1)], reference=1000.0)
    [validity] = [r for r in reasons if "validity window" in r]
    assert "too STIFF" in validity and "too SOFT" in validity


def test_the_refusal_names_the_arbiter_rather_than_asserting_a_number() -> None:
    """The window is the repo's own, and the message has to say so or it reads as a session's opinion."""
    [validity] = [r for r in refuse_unless_admissible([_pt(1)], reference=1.0)
                  if "validity window" in r]
    assert "working_stroke_strain" in validity and "MASTER-knob arbiter" in validity


def test_a_sweep_entirely_inside_the_window_is_not_blocked_by_this_check() -> None:
    """The guard must not be a blanket refusal — k_xb in [25, 100] is 5-20 nm and admissible."""
    assert _outside_working_stroke([_pt(100), _pt(50), _pt(25)]) == []
