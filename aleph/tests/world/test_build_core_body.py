"""CORE_BODY builders — the lamin meshwork and the chromatin polymer, host arithmetic only.

Each module's own ``_demo()`` is the primary check and carries the cross-checks against the source
table; these run it, because a self-check reachable only through ``python -m`` is a self-check nobody
runs.  The extra cases here are the two the PI has to be able to see failing: that a spacing from the
WRONG structure is refused rather than accepted, and that the two counts derived from one source
table still agree with the two that table reports independently.

No device is touched — the builds refuse a device-less arena, which is itself one of the assertions.
"""

from __future__ import annotations

import pytest

from aleph.world.arena import Kind, WorldArena
from aleph.world.bond import BondCount
from aleph.world.build import chromatin as chromatin_mod
from aleph.world.build import lamina as lamina_mod
from aleph.world.build.chromatin import CHROMATIN_SOURCE, SUBUNIT_UM, plan_chromatin
from aleph.world.build.envelope import nuclear_radius_um
from aleph.world.build.lamina import LAMIN_FILAMENT_UM, LAMIN_NODE_BAND, plan_lamina


def test_lamina_self_check_passes() -> None:
    """The module's own ``_demo`` is the primary check."""
    lamina_mod._demo()


def test_chromatin_self_check_passes() -> None:
    """The module's own ``_demo`` is the primary check."""
    chromatin_mod._demo()


def test_the_lamin_count_and_spacing_check_each_other() -> None:
    """Two columns of one table, taken independently: each must land where the other says.

    This is the whole argument for replacing the envelope's cortical-mesh analogy. If it ever stops
    holding, the replacement is no longer derived and the docstring claiming it is stale.
    """
    inv = plan_lamina(radius_um=nuclear_radius_um(), filament_um=LAMIN_FILAMENT_UM)
    assert LAMIN_NODE_BAND[0] <= inv["n_nodes"] <= LAMIN_NODE_BAND[1]
    assert abs(inv["areal_spacing_um_realised"] - LAMIN_FILAMENT_UM) < 1e-4


def test_the_cortical_spacing_is_refused_at_the_nuclear_radius() -> None:
    """The substitution ``envelope.ENVELOPE_MESH_PI_GAP`` makes must not pass this builder silently."""
    with pytest.raises(ValueError, match="outside the independently reported band"):
        plan_lamina(radius_um=nuclear_radius_um(), filament_um=0.100)


def test_chromatin_reports_the_shell_link_discrepancy_rather_than_hiding_it() -> None:
    """``Nc`` reproduces from its interpretation and ``Ns`` does not; both must be visible."""
    count = BondCount(basis="explicit", value=552.0, scope="HeLa MODEL value, not MCF7",
                      source_class="PI_GAP", provenance=CHROMATIN_SOURCE)
    inv = plan_chromatin(count=count, subunit_um=SUBUNIT_UM, radius_um=nuclear_radius_um())
    res = chromatin_mod._resolved(inv, count.resolve(1.0))
    assert res["n_crosslinks_derived"] == 55
    assert res["n_shell_links_implied_by_interpretation"] != inv["n_shell_links_declared"]


def test_neither_builder_has_a_host_construction_path() -> None:
    """Both refuse a device-less arena, and a refused build claims nothing."""
    arena = WorldArena(capacity={Kind.NODE: 5_000, Kind.STRAND: 4, Kind.SEGMENT: 5_000})
    count = BondCount(basis="explicit", value=552.0, scope="HeLa MODEL value, not MCF7",
                      source_class="PI_GAP", provenance=CHROMATIN_SOURCE)
    with pytest.raises(RuntimeError, match="no device allocation"):
        lamina_mod.build_lamina(arena, radius_um=nuclear_radius_um(),
                                filament_um=LAMIN_FILAMENT_UM,
                                spacing_provenance=lamina_mod.LAMIN_MESHWORK_SOURCE)
    with pytest.raises(RuntimeError, match="no device allocation"):
        chromatin_mod.build_chromatin(arena, count=count, subunit_um=SUBUNIT_UM,
                                      radius_um=nuclear_radius_um())
    assert arena.n_live(Kind.NODE) == 0
