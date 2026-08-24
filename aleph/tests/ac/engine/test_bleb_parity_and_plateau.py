"""A parity control is not a sweep point, and the softest point is not automatically a limit.

Both criteria here were declared BEFORE the runs they judge landed (2026-07-29), which is the only reason
they are worth anything: a plateau threshold chosen after seeing the ratio is not a threshold.

**Parity.** The cross-build re-run writes into the same directory as the sweep, so without help it arrives as
a second point at the reference ``k_xb`` — making the denominator ambiguous and counting the build difference
being MEASURED as a build difference that BLOCKS. Split out, it does the job it was run for: it lifts the
build-provenance refusal by measurement instead of by an operator assertion. What it must NOT do is lift more
than it covers — a pair covering two commits says nothing directly about a third, and the code has to say so.

**Plateau.** ``R`` is defined at ``k_xb -> 0``. The sweep stops at the softest point it could afford, which is
a different thing. If gamma is still falling there, the ratio BOUNDS ``R`` from above and the band comparison
is inconclusive in both directions — reading INSIDE/OUTSIDE off it is the same error as reading a mean off a
transient, the error this whole file exists to prevent one level up.
"""

from __future__ import annotations

import math

import pytest

from aleph.scripts.ac_bleb_dose_response import (
    PARITY_SIGMA,
    PLATEAU_SIGMA,
    assess_parity,
    dose_response,
    split_parity_controls,
)


def _pt(name: str, k_xb: float, mean: float, sem: float, build: str = "aaaaaaaa") -> dict:
    """A loaded record carrying only the fields these two functions read."""
    return {
        "name": name, "k_xb": k_xb, "build": build, "steps_run": 800,
        "stationarity": {"gamma_total_pn_per_um": {
            "verdict": "STATIONARY", "mean": mean, "sem": sem, "drift_over_std": 0.2}},
    }


class TestSplit:
    def test_a_repeat_of_the_same_k_xb_becomes_a_control_not_a_point(self) -> None:
        """Otherwise the reference is ambiguous — two candidates to divide by."""
        sweep, controls = split_parity_controls(
            [_pt("kxb_1000", 1000, 4.32, 0.001), _pt("kxb_1000_rebuild", 1000, 4.32, 0.001),
             _pt("kxb_10", 10, 2.88, 0.003)])
        assert [p["name"] for p in sweep] == ["kxb_1000", "kxb_10"]
        assert [p["name"] for p in controls] == ["kxb_1000_rebuild"]

    def test_the_sweep_keeps_the_point_it_always_had(self) -> None:
        """The ORIGINAL stays; the re-run is the control. Reversing it would silently restate history."""
        sweep, controls = split_parity_controls(
            [_pt("kxb_1000", 1000, 4.0, 0.01), _pt("kxb_1000_rebuild", 1000, 9.9, 0.01)])
        assert sweep[0]["name"] == "kxb_1000" and sweep[0]["stationarity"][
            "gamma_total_pn_per_um"]["mean"] == 4.0

    def test_a_sweep_with_no_repeats_is_untouched(self) -> None:
        """Backward compatibility: the three-point sweep must behave exactly as before parity existed."""
        pts = [_pt("a", 1000, 4.3, 0.001), _pt("b", 100, 3.5, 0.002), _pt("c", 10, 2.9, 0.003)]
        sweep, controls = split_parity_controls(pts)
        assert sweep == pts and controls == []


class TestParity:
    def test_agreement_is_judged_against_the_combined_sem(self) -> None:
        """Two sems add in quadrature; using one alone would call an agreeing pair a disagreement."""
        sweep, controls = split_parity_controls(
            [_pt("orig", 1000, 4.3000, 0.0010, "aaaa1111"), _pt("redo", 1000, 4.3020, 0.0010, "bbbb2222")])
        [v] = assess_parity(sweep, controls)
        assert v["agrees"] is True
        assert v["sigma"] == pytest.approx(0.0020 / math.hypot(0.001, 0.001), rel=1e-9)
        assert v["builds_covered"] == ["aaaa1111", "bbbb2222"]

    def test_a_real_disagreement_is_reported_as_one(self) -> None:
        """The control must be able to FAIL, or its agreement means nothing."""
        sweep, controls = split_parity_controls(
            [_pt("orig", 1000, 4.30, 0.001), _pt("redo", 1000, 4.40, 0.001)])
        [v] = assess_parity(sweep, controls)
        assert v["agrees"] is False and v["sigma"] > PARITY_SIGMA

    def test_a_control_with_no_mean_does_not_silently_agree(self) -> None:
        """A detector that declined to return a mean must not read as parity."""
        sweep, controls = split_parity_controls(
            [_pt("orig", 1000, 4.30, 0.001), _pt("redo", 1000, 4.30, 0.001)])
        controls[0]["stationarity"]["gamma_total_pn_per_um"]["mean"] = None
        [v] = assess_parity(sweep, controls)
        assert v["agrees"] is False and "nothing to compare" in v["why"]


class TestPlateau:
    def _block(self, soft_mean: float, soft_sem: float, next_mean: float, next_sem: float) -> dict:
        return dose_response(
            [_pt("ref", 1000, 4.3256, 0.0006), _pt("mid", 100, next_mean, next_sem),
             _pt("soft", 10, soft_mean, soft_sem)], reference=1000.0)

    def test_a_still_falling_gamma_makes_the_band_inconclusive(self) -> None:
        """The measured 168-sigma case: R bounds the limit from above and decides nothing about the band."""
        block = self._block(2.8839, 0.0031, 3.4807, 0.0017)
        assert block["softest_point_is_a_limit"] is False
        assert block["inside_band"] is None
        assert "UPPER BOUND" in block["band_verdict"]
        assert block["R"] == pytest.approx(2.8839 / 4.3256, rel=1e-12)

    def test_a_plateau_restores_the_band_verdict(self) -> None:
        """When the two softest points agree, the softest IS the limit and INSIDE/OUTSIDE is meaningful.

        1.0 / 4.3256 = 0.231, inside 0.1-0.5.  Picked against the reference rather than by eye: the first
        draft used 0.40, which is R = 0.093 and lands just BELOW the band — the test failed and was right to.
        """
        block = self._block(1.0000, 0.0100, 1.0010, 0.0100)
        assert block["softest_point_is_a_limit"] is True
        assert block["band_verdict"] == "INSIDE" and block["inside_band"] is True

    def test_a_plateau_outside_the_band_still_says_outside(self) -> None:
        """Settling is what the plateau test gates — not which side of the band the answer lands on."""
        block = self._block(3.4800, 0.0100, 3.4810, 0.0100)
        assert block["softest_point_is_a_limit"] is True and block["band_verdict"] == "OUTSIDE"

    def test_the_plateau_threshold_is_the_declared_one(self) -> None:
        """Just inside PLATEAU_SIGMA plateaus, just outside does not — no drift in the boundary."""
        sem = 0.01
        combined = math.hypot(sem, sem)
        assert self._block(1.0, sem, 1.0 + (PLATEAU_SIGMA - 0.1) * combined, sem)[
            "softest_point_is_a_limit"] is True
        assert self._block(1.0, sem, 1.0 + (PLATEAU_SIGMA + 0.1) * combined, sem)[
            "softest_point_is_a_limit"] is False

    def test_the_note_states_the_measured_step_so_the_verdict_is_auditable(self) -> None:
        """A boolean alone cannot be checked; the sigma that produced it must travel with it."""
        note = self._block(2.8839, 0.0031, 3.4807, 0.0017)["plateau_note"]
        assert "sigma" in note and "0.5968" in note
