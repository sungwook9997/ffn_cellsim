"""The cytosol module refuses, and the refusal is checkable rather than only readable.

Each module's own ``_demo()`` is the primary check (the ``world`` convention); this file is the thin
pytest collector for it, mirroring ``test_build_core_body.py``.

What it adds beyond ``_demo()`` is the one property Lead asked for explicitly: while
``arena.Kind`` has no ``GRID_CELL``, a caller must get a TYPED, self-describing refusal and never an
``AttributeError``.  The difference is not cosmetic — an ``AttributeError`` reads as a defect in this
file, and the actual state is a decision pending in another one.  When the ``Kind`` lands these tests
change meaning without changing text: the refusal branch stops applying and ``_demo`` skips it.
"""

from __future__ import annotations

from inspect import Parameter, signature

import pytest

from aleph.world import build as build_pkg  # noqa: F401  (package import must not drag in a device)
from aleph.world.arena import Kind
from aleph.world.build import cytosol as cytosol_mod


def test_module_demo_passes() -> None:
    """The module's own ``_demo`` is the primary check; a failure there is a failure here."""
    cytosol_mod._demo()


@pytest.mark.skipif(hasattr(Kind, "GRID_CELL"), reason="Kind.GRID_CELL landed; the gap is closed")
def test_missing_grid_kind_is_a_typed_refusal_not_an_attribute_error() -> None:
    """A pending decision must be distinguishable from a broken module, by type and by message."""
    with pytest.raises(cytosol_mod.ArenaKindGapError) as excinfo:
        cytosol_mod.grid_cell_kind()
    err = excinfo.value
    assert not isinstance(err, AttributeError)
    assert err.kind_name == "GRID_CELL"
    assert err.population == "cytosol"
    # The refusal carries its own argument and where the decision is recorded, so a reader who finds
    # it in a traceback does not have to go looking for the reason.
    assert "arena.py:77-78" in err.why
    assert "CYTOSOL_ARENA_REPRESENTATION_2026-08-20.md" in err.why


def test_the_pi_resolution_is_recorded_but_is_not_a_default() -> None:
    """PI set ``dx_um`` 2026-08-20; a decision is recorded, and recording it must not create a default.

    The distinction is the whole point of the axis table. A value in ``CYTOSOL_AXES`` says *someone
    decided this, here is who and why*; a value in the SIGNATURE would let a run inherit a resolution
    it never stated, and ``dx`` moves the fluid cost as dx^-5. So the axis carries 0.25 and
    :func:`build_cytosol` still refuses to be called without one.
    """
    axis = cytosol_mod.CYTOSOL_AXES["dx_um"]
    assert axis["value"] == cytosol_mod.PI_DX_UM
    assert "PI-DECISION" in str(axis["provenance"])
    assert "test point" in str(axis["provenance"]), "a test point must not read as a derivation"
    assert signature(cytosol_mod.build_cytosol).parameters["dx_um"].default is Parameter.empty

    # The axes with no PI point still carry none, and every axis carries its provenance text.
    assert cytosol_mod.CYTOSOL_AXES["poroelastic_diffusion_um2_s"]["value"] is None
    assert cytosol_mod.CYTOSOL_AXES["dt_phys_s"]["value"] is None
    for name, a in cytosol_mod.CYTOSOL_AXES.items():
        assert a["source"], f"{name} has no provenance text"


def test_both_ends_of_the_sourced_band_refuse() -> None:
    """The band is enforced, not only described — and the PI point sits inside it."""
    common = dict(radius_um=7.5, poroelastic_diffusion_um2_s=50.0, dt_phys_s=0.05)
    with pytest.raises(ValueError, match="continuum floor"):
        cytosol_mod.cytosol_counts(dx_um=0.005, **common)       # below the pore size
    with pytest.raises(ValueError, match="diffusion length"):
        cytosol_mod.cytosol_counts(dx_um=3.0, **common)         # coarser than one step's transport
    counts = cytosol_mod.cytosol_counts(dx_um=cytosol_mod.PI_DX_UM, **common)
    assert counts["n_cells"] == 63 ** 3
