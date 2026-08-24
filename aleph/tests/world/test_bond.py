"""Bond families — the count that cannot be a bare number, and the pair that is derived.

The property under test is not arithmetic. It is that a family CANNOT be constructed without answering
how many, against what, for which cell, and on whose authority — because the incumbent's ERM count
became a mesh artefact precisely by having a fallback for that question.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.units.provenance import PROVENANCE
from aleph.world.arena import Kind, WorldArena
from aleph.world.bond import BondCount, SourceClass, build_bond_family
from aleph.world.strand import build_strand
from aleph.world import bond as bond_mod


def _count(**over) -> BondCount:
    base = dict(basis="areal", value=10.0, scope="MCF7 interphase adherent 37C",
                source_class=SourceClass.PI_GAP, provenance="no Contract-Graph datum")
    return BondCount(**{**base, **over})


@pytest.fixture()
def world() -> tuple[WorldArena, np.ndarray]:
    """An arena with two parallel strands and the candidate pair list between them."""
    arena = WorldArena(capacity={Kind.NODE: 500, Kind.SEGMENT: 500, Kind.ANGLE3: 500, Kind.BOND: 500})
    a = build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0), contour_um=1.0, seg_um=0.1)
    b = build_strand(arena, "sf_arc", start=(0, 0.05, 0), direction=(1, 0, 0),
                     contour_um=1.0, seg_um=0.1)
    pairs = np.stack([a.nodes.lo + np.arange(a.n_nodes), b.nodes.lo + np.arange(b.n_nodes)], axis=1)
    return arena, pairs


def test_self_check_passes() -> None:
    bond_mod._demo()


def test_source_class_has_not_drifted_from_the_engine_vocabulary() -> None:
    """One vocabulary on purpose — a second spelling lets a typo pass silently."""
    assert set(PROVENANCE) <= {c.value for c in SourceClass}
    assert "UNRATIFIED_PROXY" in {c.value for c in SourceClass}, "the scope-mismatch class must exist"


def test_scope_is_required() -> None:
    """A sourced value out of scope is a proxy: Pi_0 = 40 Pa is a real measurement, of HeLa."""
    with pytest.raises(ValueError, match="scope is required"):
        _count(scope="   ")


def test_provenance_is_required() -> None:
    """A count with no stated origin is how a discretisation artefact acquires a physiological label."""
    with pytest.raises(ValueError, match="provenance is required"):
        _count(provenance="")


@pytest.mark.parametrize("bad", ["per_angstrom", "", "AREAL"])
def test_basis_must_be_known(bad: str) -> None:
    with pytest.raises(ValueError, match="basis must be one of"):
        _count(basis=bad)


@pytest.mark.parametrize("bad", [-1.0, float("nan"), float("inf")])
def test_value_must_be_finite_and_nonnegative(bad: float) -> None:
    with pytest.raises(ValueError, match="finite and nonnegative"):
        _count(value=bad)


def test_density_resolves_against_a_measured_support() -> None:
    """The support is the measured mesh area, never the analytic sphere — that is what keeps an areal
    density invariant when the continuum resolution changes."""
    c = _count(basis="areal", value=58.0)
    assert c.support_unit == "um^2"
    assert c.resolve(10.0) == 580
    assert c.resolve(0.5) == 29


def test_zero_bonds_is_an_answer_not_an_error() -> None:
    """A density low enough to give no bonds over the support is physical, unlike a count that was
    never asked for."""
    assert _count(value=1e-9).resolve(1.0) == 0


@pytest.mark.parametrize("bad", [0.0, -2.0, float("nan")])
def test_support_must_be_positive(bad: float) -> None:
    with pytest.raises(ValueError, match="finite and positive"):
        _count().resolve(bad)


def test_explicit_basis_takes_no_support() -> None:
    c = _count(basis="explicit", value=7.0)
    assert c.resolve(1.0) == 7
    with pytest.raises(ValueError, match="takes no support"):
        c.resolve(4.0)


def test_component_pair_is_derived_not_declared(world) -> None:
    """A bond stores two global node ids; which populations it joins is read off the ID ranges."""
    arena, pairs = world
    fam = build_bond_family(arena, "sf_cortex_transient", chemistry_card="transient_actin",
                            count=_count(basis="explicit", value=6.0), support=1.0,
                            pairs=pairs, rest_um=0.05, stiffness_pn_per_um=250.0)
    assert fam.component_pairs(arena) == {("cortex", "sf_arc"): 6}


def test_internal_bond_needs_no_separate_declaration() -> None:
    """A bond inside one population reports that population twice — no INTERNAL/INTER split needed."""
    arena = WorldArena(capacity={Kind.NODE: 200, Kind.SEGMENT: 200, Kind.ANGLE3: 200, Kind.BOND: 200})
    a = build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0), contour_um=1.0, seg_um=0.1)
    b = build_strand(arena, "cortex", start=(0, 0.05, 0), direction=(1, 0, 0), contour_um=1.0, seg_um=0.1)
    pairs = np.stack([a.nodes.lo + np.arange(a.n_nodes), b.nodes.lo + np.arange(b.n_nodes)], axis=1)
    fam = build_bond_family(arena, "alpha_actinin", chemistry_card="alpha_actinin_ferrer2008",
                            count=_count(basis="explicit", value=5.0), support=1.0,
                            pairs=pairs, rest_um=0.05, stiffness_pn_per_um=4.6e5)
    assert fam.component_pairs(arena) == {("cortex", "cortex"): 5}


def test_provenance_row_keeps_count_and_scope_together(world) -> None:
    """Together on purpose: a count separated from its scope is how "235/um^2" came to look like a
    measurement when it was the subdiv-7 vertex density."""
    arena, pairs = world
    fam = build_bond_family(arena, "erm", chemistry_card="ezrin",
                            count=_count(basis="areal", value=58.0,
                                         source_class=SourceClass.UNRATIFIED_PROXY,
                                         provenance="subdiv-6 vertex density, NOT physiological"),
                            support=0.1, pairs=pairs, rest_um=0.02, stiffness_pn_per_um=4.6e3)
    row = fam.provenance_row()
    for key in ("scope", "provenance", "source_class", "basis", "support", "support_unit", "n_bonds"):
        assert row[key] not in (None, ""), key
    assert row["source_class"] == "UNRATIFIED_PROXY"


def test_too_few_candidates_refuses_rather_than_truncating(world) -> None:
    """Supply more candidates, never lower the density: the density is a physiological axis and the
    candidate list is a build detail."""
    arena, pairs = world
    with pytest.raises(ValueError, match="Supply more candidates"):
        build_bond_family(arena, "erm", chemistry_card="ezrin",
                          count=_count(basis="explicit", value=float(pairs.shape[0] + 1)),
                          support=1.0, pairs=pairs, rest_um=0.02, stiffness_pn_per_um=1.0)


def test_bond_to_an_unclaimed_node_is_refused(world) -> None:
    arena, _ = world
    beyond = arena.n_live(Kind.NODE) + 3
    with pytest.raises(AssertionError, match="connection to nothing"):
        build_bond_family(arena, "ghost", chemistry_card="none",
                          count=_count(basis="explicit", value=1.0), support=1.0,
                          pairs=np.array([[0, beyond]]), rest_um=0.0, stiffness_pn_per_um=1.0)


def test_self_loop_is_refused(world) -> None:
    arena, _ = world
    with pytest.raises(ValueError, match="no direction"):
        build_bond_family(arena, "loop", chemistry_card="none",
                          count=_count(basis="explicit", value=1.0), support=1.0,
                          pairs=np.array([[3, 3]]), rest_um=0.0, stiffness_pn_per_um=1.0)


@pytest.mark.parametrize("rest,stiff", [(-0.1, 1.0), (float("nan"), 1.0), (0.05, -1.0)])
def test_negative_or_non_finite_mechanics_are_refused(world, rest: float, stiff: float) -> None:
    arena, pairs = world
    with pytest.raises(ValueError, match="finite and nonnegative"):
        build_bond_family(arena, "bad", chemistry_card="none",
                          count=_count(basis="explicit", value=2.0), support=1.0,
                          pairs=pairs, rest_um=rest, stiffness_pn_per_um=stiff)


def test_family_claims_a_bond_range_and_the_arena_stays_partitioned(world) -> None:
    arena, pairs = world
    fam = build_bond_family(arena, "erm", chemistry_card="ezrin",
                            count=_count(basis="explicit", value=4.0), support=1.0,
                            pairs=pairs, rest_um=0.02, stiffness_pn_per_um=4.6e3)
    assert fam.bonds.count == 4 and arena.n_live(Kind.BOND) == 4
    arena.assert_partitioned()
    assert arena.census()["populations"]["erm"]["bond"] == 4
