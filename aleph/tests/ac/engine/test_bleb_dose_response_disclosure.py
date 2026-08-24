"""What the STATIONARY label hides must reach the reader.

The dose-response analyser refuses to compute R unless every point settled.  That refusal is only as good
as the thing it reads, and it reads exactly one observable's verdict (``gamma_total``, correctly — R is a
``gamma_total`` ratio).  On the 2026-07-29 sweep that left two true things invisible:

* ``k_xb=10`` passed at 0.98 sigma against a 1.0 threshold — decided by 2%, printed identically to a point
  that passed at 0.13;
* its ``gamma_source`` did NOT settle, and the detector DECLINED TO RETURN A MEAN for it.

Neither should block R.  Both must be printed.  A refusal the reader never sees is the same defect the
analyser exists to prevent, one level down — which is why this is tested rather than left to habit.
"""

from __future__ import annotations

import pytest

from aleph.scripts.ac_bleb_dose_response import MARGINAL_DRIFT, _disclose_marginal_and_drifting


def _point(name: str, *, drift: float, source_verdict: str = "STATIONARY",
           source_drift: float | None = 0.2) -> dict:
    """One loaded sweep point, carrying only the fields the disclosure reads."""
    return {
        "name": name,
        "stationarity": {
            "gamma_total_pn_per_um": {"verdict": "STATIONARY", "mean": 3.0, "sem": 0.01,
                                      "drift_over_std": drift},
            "gamma_source_pn_per_um": {"verdict": source_verdict,
                                       "mean": None if source_verdict != "STATIONARY" else 0.02,
                                       "drift_over_std": source_drift},
        },
    }


def test_a_comfortable_pass_says_nothing(capsys: pytest.CaptureFixture[str]) -> None:
    """Silence is the correct output when there is nothing the label hides."""
    _disclose_marginal_and_drifting([_point("kxb_1000", drift=0.13)])
    assert capsys.readouterr().out == ""


def test_a_marginal_pass_is_named_with_its_margin(capsys: pytest.CaptureFixture[str]) -> None:
    """0.98 against 1.0 is a pass, and the reader must be told by how little."""
    _disclose_marginal_and_drifting([_point("kxb_10", drift=0.98)])
    out = capsys.readouterr().out
    assert "kxb_10" in out and "0.98" in out and "marginal" in out


def test_the_threshold_itself_is_the_boundary(capsys: pytest.CaptureFixture[str]) -> None:
    """Exactly at MARGINAL_DRIFT discloses; just below stays silent — no off-by-one in the guard."""
    _disclose_marginal_and_drifting([_point("at", drift=MARGINAL_DRIFT)])
    assert "marginal" in capsys.readouterr().out
    _disclose_marginal_and_drifting([_point("below", drift=MARGINAL_DRIFT - 0.01)])
    assert capsys.readouterr().out == ""


def test_a_drifting_sibling_is_disclosed_but_does_not_read_as_fatal(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """It must be visible AND must say why it is not a blocking reason, or a reader will treat it as one."""
    _disclose_marginal_and_drifting(
        [_point("kxb_10", drift=0.5, source_verdict="DRIFTING", source_drift=1.18)]
    )
    out = capsys.readouterr().out
    assert "gamma_source_pn_per_um" in out and "did NOT settle" in out
    assert "does not block" in out


def test_a_sibling_with_no_recorded_drift_still_reports(capsys: pytest.CaptureFixture[str]) -> None:
    """A missing number must not silence the disclosure — that is how the failure would hide."""
    _disclose_marginal_and_drifting(
        [_point("kxb_x", drift=0.5, source_verdict="TOO_SHORT", source_drift=None)]
    )
    assert "no drift recorded" in capsys.readouterr().out


def test_every_point_is_examined_not_just_the_first(capsys: pytest.CaptureFixture[str]) -> None:
    """The sweep is a list; disclosing only the first point would hide the softest, which is the one at risk."""
    _disclose_marginal_and_drifting(
        [_point("kxb_1000", drift=0.13), _point("kxb_100", drift=0.19), _point("kxb_10", drift=0.98)]
    )
    out = capsys.readouterr().out
    assert "kxb_10" in out and "kxb_1000" not in out
