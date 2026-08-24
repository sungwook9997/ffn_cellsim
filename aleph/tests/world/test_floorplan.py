"""Floorplan — the page must show what is wired, and must not flatter it.

The properties under test are the two the view exists for: that a count nobody sourced does not render
like one that was sourced, and that the page is self-contained. Everything else is layout.
"""

from __future__ import annotations

import numpy as np
import pytest

from aleph.world.arena import Kind, WorldArena
from aleph.world.bond import BondCount, SourceClass, build_bond_family
from aleph.world.floorplan import render_floorplan, write_floorplan
from aleph.world.strand import build_strand
from aleph.world import floorplan as fp_mod


def _built() -> tuple[WorldArena, list]:
    arena = WorldArena(capacity={Kind.NODE: 5_000, Kind.SEGMENT: 5_000,
                                 Kind.ANGLE3: 5_000, Kind.BOND: 5_000})
    cortex = [build_strand(arena, "cortex", start=(0, 0.05 * i, 0), direction=(1, 0, 0),
                           contour_um=1.0, seg_um=0.1) for i in range(6)]
    memb = build_strand(arena, "membrane", start=(0, 0, 1.0), direction=(1, 0, 0),
                        contour_um=1.0, seg_um=0.1)
    cn = np.concatenate([s.nodes.lo + np.arange(s.n_nodes) for s in cortex])
    mn = memb.nodes.lo + np.arange(memb.n_nodes)
    fams = [
        build_bond_family(arena, "alpha_actinin", chemistry_card="alpha_actinin_ferrer2008",
                          count=BondCount(basis="explicit", value=12.0, scope="MCF7 interphase 37C",
                                          source_class=SourceClass.SOURCED,
                                          provenance="Ferrer 2008 PNAS AFM; PI-approved 2026-06-30"),
                          support=1.0, pairs=np.stack([cn[:-1], cn[1:]], axis=1),
                          rest_um=0.05, stiffness_pn_per_um=4.6e5),
        build_bond_family(arena, "erm", chemistry_card="ezrin",
                          count=BondCount(basis="areal", value=58.0, scope="none — subdiv-6 mesh",
                                          source_class=SourceClass.UNRATIFIED_PROXY,
                                          provenance="subdiv-6 vertex density; NOT physiological"),
                          support=0.1, pairs=np.stack([mn, cn[:mn.size]], axis=1),
                          rest_um=0.02, stiffness_pn_per_um=4.6e3),
    ]
    return arena, fams


def test_self_check_passes() -> None:
    fp_mod._demo()


def test_empty_arena_renders_rather_than_raising() -> None:
    """"Nothing is wired yet" is the state this view exists to show at the start of a build."""
    page = render_floorplan(WorldArena(capacity={Kind.NODE: 100}))
    assert "nothing is wired" in page
    assert page.startswith("<!doctype html>")


def test_page_is_self_contained() -> None:
    """No external asset of any kind — the page must render with no network at all."""
    arena, fams = _built()
    page = render_floorplan(arena, fams)
    for forbidden in ("http://", "https://", "<script src", "<link ", "@import"):
        assert forbidden not in page, forbidden


def test_unsourced_counts_are_flagged_and_counted() -> None:
    """The whole point: a count nobody sourced must not look like one that was sourced."""
    arena, fams = _built()
    page = render_floorplan(arena, fams)
    assert "1 of 2 families rest on a count that is not a sourced measurement" in page
    assert "UNRATIFIED_PROXY" in page and "NOT physiological" in page
    assert "Ferrer 2008" in page, "the sourced family's authority must be shown too"


def test_all_sourced_says_so_instead_of_staying_silent() -> None:
    """Silence would be indistinguishable from "the panel is broken"."""
    arena = WorldArena(capacity={Kind.NODE: 200, Kind.SEGMENT: 200, Kind.ANGLE3: 200, Kind.BOND: 200})
    s = build_strand(arena, "cortex", start=(0, 0, 0), direction=(1, 0, 0), contour_um=1.0, seg_um=0.1)
    n = s.nodes.lo + np.arange(s.n_nodes)
    fam = build_bond_family(arena, "alpha_actinin", chemistry_card="a",
                            count=BondCount(basis="explicit", value=3.0, scope="MCF7",
                                            source_class=SourceClass.SOURCED, provenance="Ferrer 2008"),
                            support=1.0, pairs=np.stack([n[:-1], n[1:]], axis=1),
                            rest_um=0.05, stiffness_pn_per_um=1.0)
    page = render_floorplan(arena, [fam])
    assert "Every family's count is SOURCED or DERIVED" in page


def test_totals_match_the_census_exactly() -> None:
    """The page cannot disagree with the arena, because it reads the arena's own census."""
    arena, fams = _built()
    page = render_floorplan(arena, fams)
    assert f"{arena.n_live(Kind.NODE):,} nodes live" in page
    assert f"{arena.n_live(Kind.BOND):,}" in page


def test_derived_joins_are_shown_and_were_never_declared() -> None:
    """cortex-membrane appears on the page although no edge anywhere declares it."""
    arena, fams = _built()
    page = render_floorplan(arena, fams)
    assert "cortex&ndash;membrane" in page or "cortex–membrane" in page


def test_scope_and_provenance_travel_with_the_count() -> None:
    """A count separated from its scope is how a mesh number came to look like a measurement."""
    arena, fams = _built()
    page = render_floorplan(arena, fams)
    for f in fams:
        assert f.count.scope.split(" —")[0] in page
        assert f.count.provenance.split(";")[0] in page


def test_write_floorplan_round_trips(tmp_path) -> None:
    arena, fams = _built()
    out = write_floorplan(arena, str(tmp_path / "fp.html"), fams, title="round trip")
    text = open(out, encoding="utf-8").read()
    assert "round trip" in text and text.rstrip().endswith("</html>")


def test_title_and_names_are_escaped() -> None:
    """A population name is caller data and must not be able to inject markup."""
    arena = WorldArena(capacity={Kind.NODE: 100, Kind.SEGMENT: 100, Kind.ANGLE3: 100})
    build_strand(arena, "<script>x</script>", start=(0, 0, 0), direction=(1, 0, 0),
                 contour_um=1.0, seg_um=0.1)
    page = render_floorplan(arena, title="<b>t</b>")
    assert "<script>x</script>" not in page and "&lt;script&gt;" in page
    assert "<title><b>t</b></title>" not in page
