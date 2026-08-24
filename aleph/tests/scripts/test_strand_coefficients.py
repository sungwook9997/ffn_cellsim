"""Coefficients must be derived from sourced constants, or refused — never defaulted.

⚠ **WHY THIS EXISTS.** Eleven populations stand in the arena and two carry force. The seven pairs of
numbers that would change that are what this file guards, and the guard is not "are they right" — it
is **can a plausible number reach a population without a derivation behind it**. A microtubule handed
actin's rigidity would run, look fine, and be a cell made of the wrong material.
"""

from __future__ import annotations

import math

import pytest

from aleph.laws import units as U
from aleph.scripts.strand_coefficients import (
    D_IF_UM,
    EA_ACTIN_PN,
    E_IF_PN_UM2,
    K_CHROMATIN_PN_PER_UM,
    K_CORTEX_PN_PER_UM,
    LP_IF_BAND_UM,
    R_MT_LUMEN_UM,
    R_MT_OUTER_UM,
    derive,
    kappa_from_persistence_length,
    kappa_from_rod,
)


def _by():
    return {c.population: c for c in derive()}


def test_a_microtubule_never_gets_actins_rigidity() -> None:
    """The single failure that would be invisible in a run and fatal to every conclusion from it."""
    by = _by()
    assert by["microtubule"].kappa_pn_um2 == U.KAPPA_MT
    assert by["microtubule"].kappa_pn_um2 > 100.0 * by["stress_fiber"].kappa_pn_um2


def test_the_mt_lumen_radius_is_licensed_by_a_cross_check_not_by_choice() -> None:
    """⚠ The lumen radius is only a derivation while E stays inside the literature range.

    It was declared 2026-08-24 to unblock EA. What stops it being a free knob is that the ALREADY
    SOURCED bending rigidity over-determines it: `E = KAPPA_MT / I` must land in 1-2 GPa. A wall
    thickness chosen to make something else pass would push E out of that range and fail here.
    """
    i_mt = math.pi * (R_MT_OUTER_UM**4 - R_MT_LUMEN_UM**4) / 4.0
    e_gpa = (U.KAPPA_MT / i_mt) / 1e9
    assert 1.0 <= e_gpa <= 2.0, e_gpa
    assert R_MT_LUMEN_UM < R_MT_OUTER_UM


def test_the_two_independent_routes_to_if_rigidity_agree() -> None:
    """IF is over-determined: an elastic rod from E and d, and a persistence-length band.

    They come from different papers (Kreplak 2005 / Guo 2013 for E; Mücke 2004 for d; the band from a
    third). Agreement is what makes the value DERIVED. If they ever stop agreeing, neither may be
    used, and the test says so rather than quietly averaging them.
    """
    rod = kappa_from_rod(E_IF_PN_UM2, D_IF_UM)
    lo, hi = (kappa_from_persistence_length(x) for x in LP_IF_BAND_UM)
    assert lo <= rod <= hi, (lo, rod, hi)


def test_one_material_has_one_stiffness_across_populations() -> None:
    """⚠ Actin is actin. The cortex ran at 88,000 while stress fibres would take 880,000.

    Ratified 2026-08-24. Leaving them different lets the LABEL rather than the MATERIAL set the
    physics, which is the same class the charter's "one filament, one component" rule forbids.
    """
    by = _by()
    actin = {"cortex", "stress_fiber", "sf_arc", "lamellipodium", "filopodium"}
    stiffnesses = {by[p].k_axial_pn_per_um for p in actin}
    assert len(stiffnesses) == 1, {p: by[p].k_axial_pn_per_um for p in actin}
    assert by["cortex"].k_axial_pn_per_um == K_CORTEX_PN_PER_UM == EA_ACTIN_PN / 0.05
    kappas = {by[p].kappa_pn_um2 for p in actin}
    assert len(kappas) == 1, kappas


def test_chromatin_carries_a_sourced_spring_and_no_bending_triple() -> None:
    """⚠ Both halves were wrong in the first draft and both are structural, not numerical.

    `kp` = 1.6 nN/µm is SOURCED (Stephens 2017 Table 1) and is already a spring constant, so it is
    NOT divided by a step. And `build/chromatin` claims no ANGLE3 range at all — the source models
    chromatin as extensible springs of ZERO bending modulus — so a bending term has nothing to launch
    over and adding one would add physics the source does not have.
    """
    by = _by()
    assert by["chromatin"].k_axial_pn_per_um == K_CHROMATIN_PN_PER_UM == 1.6e3
    assert by["chromatin"].kappa_pn_um2 is None
    assert by["chromatin"].flag() == "--bind chromatin=1600"   # no kappa: nothing to launch over
    assert by["chromatin"].grade == "SOURCED"


def test_the_cortex_does_not_go_through_bind() -> None:
    """⚠ `build_all` stands the cortex, so `--bind cortex=...` is REFUSED by the driver at launch.

    The first draft of `flag()` emitted it and the run would have died on
    "refused: --bind names 'cortex', which PHASE 1 did not build".
    """
    flag = _by()["cortex"].flag()
    assert flag is not None and flag.startswith("--cortex-k-axial"), flag
    assert "--bind" not in flag


def test_every_coefficient_is_positive_and_finite() -> None:
    for c in derive():
        for value in (c.kappa_pn_um2, c.k_axial_pn_per_um):
            if value is not None:
                assert math.isfinite(value) and value > 0.0, (c.population, value)


def test_a_non_positive_step_is_refused() -> None:
    """The step divides EA and cubes into bending; zero must raise, not produce infinity."""
    for bad in (0.0, -0.05):
        with pytest.raises(ValueError):
            derive(bad)


def test_axial_scales_with_the_step_and_bending_does_not() -> None:
    """κ is a MATERIAL property; k_axial is EA/L. Getting this backwards breaks grid invariance."""
    a, b = _by(), {c.population: c for c in derive(0.10)}
    assert a["stress_fiber"].kappa_pn_um2 == b["stress_fiber"].kappa_pn_um2
    assert a["stress_fiber"].k_axial_pn_per_um == pytest.approx(
        2.0 * b["stress_fiber"].k_axial_pn_per_um)
    # ⚠ Chromatin's spring is NOT a step-divided quantity, so it must not move with the step.
    assert a["chromatin"].k_axial_pn_per_um == b["chromatin"].k_axial_pn_per_um


def test_the_bundle_understatement_travels_with_the_value() -> None:
    """Filopodia and stress fibres are BUNDLES built as single chains. The value understates them.

    Cards F1-F4 hold the fascin constants that would fix it and none is bound. A reader taking these
    numbers must see that, so it is pinned to the basis string rather than left in a commit message.
    """
    by = _by()
    assert "bundle" in by["filopodium"].basis.lower()
    assert "F1-F4" in by["filopodium"].basis
    assert "BUNDLE" in by["stress_fiber"].basis or "bundle" in by["stress_fiber"].basis
