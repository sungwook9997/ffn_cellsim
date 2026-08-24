"""Cell view — every element reaches the payload, and the payload decodes back to what went in.

The property under test is the one the page claims on its own face: nothing is sampled away. A viewer
that silently dropped elements would be indistinguishable from a smaller cell, so the encode/decode
round-trip and the per-layer counts are what make the claim checkable rather than asserted.
"""

from __future__ import annotations

import base64
import json

import numpy as np
import pytest

from aleph.world.arena import Kind, WorldArena
from aleph.world.bond import BondCount, SourceClass, build_bond_family
from aleph.world.strand import build_strand
from aleph.world.surface import build_surface, icosphere
from aleph.world.viewer import render_cell_view, write_cell_view
from aleph.world import viewer as vw


def _scene(page: str) -> dict:
    return json.loads(page.split('id="scene">')[1].split("</script>")[0])


def _built():
    arena = WorldArena(capacity={Kind.NODE: 9_000, Kind.SEGMENT: 9_000, Kind.ANGLE3: 9_000,
                                 Kind.ANGLE4: 9_000, Kind.FACE: 9_000, Kind.BOND: 9_000})
    strands = [build_strand(arena, "cortex", start=(-1.5, 0.1 * i, 0), direction=(1, 0, 0),
                            contour_um=3.0, seg_um=0.1) for i in range(8)]
    v, f = icosphere(2)
    surfaces = [build_surface(arena, "membrane", vertices=v, faces=f, radius_um=2.5)]
    cn = np.concatenate([s.nodes.lo + np.arange(s.n_nodes) for s in strands])
    fams = [build_bond_family(
        arena, "alpha_actinin", chemistry_card="alpha_actinin_ferrer2008",
        count=BondCount(basis="explicit", value=60.0, scope="MCF7 interphase 37C",
                        source_class=SourceClass.SOURCED, provenance="Ferrer 2008 PNAS"),
        support=1.0, pairs=np.stack([cn[:-31], cn[31:]], axis=1),
        rest_um=0.1, stiffness_pn_per_um=4.6e5)]
    return arena, strands, surfaces, fams


def test_self_check_passes() -> None:
    vw._demo()


def test_every_element_reaches_the_payload() -> None:
    """N-1 segments per strand, unique edges per surface, every bond — none dropped."""
    arena, strands, surfaces, fams = _built()
    s = _scene(render_cell_view(arena, strands, surfaces, fams))
    by = {(l["name"], l["kind"]): l["n"] for l in s["layers"]}
    assert by[("cortex", "strand")] == sum(st.segments.count for st in strands)
    assert by[("membrane", "surface")] == 3 * surfaces[0].faces.count // 2
    assert by[("alpha_actinin", "bond")] == fams[0].n_bonds


def test_positions_round_trip_through_the_encoding() -> None:
    """float32 base64 is a transport choice, not a filter: what decodes must be what was built."""
    arena, strands, surfaces, fams = _built()
    s = _scene(render_cell_view(arena, strands, surfaces, fams))
    back = np.frombuffer(base64.b64decode(s["pos"]), np.float32).reshape(-1, 3)
    assert back.shape[0] == arena.n_live(Kind.NODE) == s["n_pos"]
    for st in strands:
        assert np.allclose(back[st.nodes.lo:st.nodes.hi], st.position, atol=1e-5)
    assert np.allclose(back[surfaces[0].nodes.lo:surfaces[0].nodes.hi], surfaces[0].position, atol=1e-5)


def test_indices_round_trip_and_stay_in_range() -> None:
    arena, strands, surfaces, fams = _built()
    s = _scene(render_cell_view(arena, strands, surfaces, fams))
    for layer in s["layers"]:
        idx = np.frombuffer(base64.b64decode(layer["idx"]), np.uint32).reshape(-1, 2)
        assert idx.shape[0] == layer["n"]
        assert idx.max() < s["n_pos"], f"{layer['name']} indexes past the position table"


def test_page_is_self_contained() -> None:
    arena, strands, surfaces, fams = _built()
    page = render_cell_view(arena, strands, surfaces, fams)
    for forbidden in ("http://", "https://", "<script src", "<link ", "@import"):
        assert forbidden not in page, forbidden


def test_source_class_travels_onto_the_layer_chip() -> None:
    """The provenance follows the geometry into the view — a PI_GAP family must say so on screen."""
    arena, strands, surfaces, _ = _built()
    cn = np.concatenate([s.nodes.lo + np.arange(s.n_nodes) for s in strands])
    gap = build_bond_family(
        arena, "nmii_crossbridge", chemistry_card="nmii_catch_slip_kovacs2007",
        count=BondCount(basis="explicit", value=10.0, scope="MCF7",
                        source_class=SourceClass.PI_GAP, provenance="cortical_minifilament_density"),
        support=1.0, pairs=np.stack([cn[:-5], cn[5:]], axis=1), rest_um=0.0, stiffness_pn_per_um=1.0)
    page = render_cell_view(arena, strands, surfaces, [gap])
    assert "PI_GAP" in page


def test_empty_arena_does_not_divide_by_a_zero_extent() -> None:
    page = render_cell_view(WorldArena(capacity={Kind.NODE: 10}))
    s = _scene(page)
    assert s["extent"] > 0.0 and s["layers"] == []


def test_extent_is_read_from_the_geometry() -> None:
    arena, strands, surfaces, fams = _built()
    s = _scene(render_cell_view(arena, strands, surfaces, fams))
    assert s["extent"] == pytest.approx(2.5, abs=0.2), "the membrane radius dominates this scene"


def test_fragment_mode_omits_the_document_skeleton() -> None:
    arena, strands, surfaces, fams = _built()
    frag = render_cell_view(arena, strands, surfaces, fams, standalone=False)
    assert not frag.startswith("<!doctype") and "<body>" not in frag and "</html>" not in frag
    assert "id=\"scene\"" in frag, "the fragment is the same page, unwrapped"


def test_note_is_shown_and_escaped() -> None:
    arena, strands, surfaces, fams = _built()
    page = render_cell_view(arena, strands, surfaces, fams, note="5.7% of native <b>x</b>")
    assert "5.7% of native" in page and "<b>x</b>" not in page


def test_write_cell_view_round_trips(tmp_path) -> None:
    arena, strands, surfaces, fams = _built()
    out = write_cell_view(arena, str(tmp_path / "c.html"), strands, surfaces, fams, title="rt")
    assert "rt" in open(out, encoding="utf-8").read()
