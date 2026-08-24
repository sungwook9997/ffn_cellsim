"""Controls for the layer restructure. Written before the move, because the rename was not.

The rename ran first, missed three positions the same day, and cost twenty-two red tests plus a
`Makefile` and a git hook that pointed at paths which no longer existed. Every one of those is a
case here, so this move pays for that lesson once rather than twice.
"""

from __future__ import annotations

import pytest

from aleph.scripts.restructure_to_layers import MOVES, _pairs, rewrite


@pytest.mark.parametrize(
    ("before", "after"),
    [
        # dotted — the ordinary import
        ("from aleph.ac.engine import contracts", "from aleph.engine import contracts"),
        ("import aleph.ff.network_warp", "import aleph.laws.network_warp"),
        ("from aleph.ac.motor.hand import Hand", "from aleph.components.motor.hand import Hand"),
        # the incumbent relocates rather than leaving
        ("from aleph.ac.cell.assemble import build", "from aleph.components.incumbent.assemble import build"),
        # slashed — a path in prose, a Makefile target, a hook
        ("aleph/ac/engine/contracts.py", "aleph/engine/contracts.py"),
        ("$(PY) aleph/ff/units.py", "$(PY) aleph/laws/units.py"),
        # component-wise quoted — the position the rename missed
        ('ROOT / "aleph" / "ac" / "cell"', 'ROOT / "aleph" / "components" / "incumbent"'),
        ('REPO_ROOT / "aleph" / "ff"', 'REPO_ROOT / "aleph" / "laws"'),
    ],
)
def test_every_spelling_moves(before: str, after: str) -> None:
    rewritten, count = rewrite(before)
    assert rewritten == after
    assert count == 1


def test_a_longer_source_wins_over_its_prefix() -> None:
    """`ac.cell` and `ac` are both sources; matching `ac` first would produce `aleph.cell`."""
    sources = [old for old, _ in _pairs()]
    assert sources == sorted(sources, key=len, reverse=True)

    rewritten, _ = rewrite("aleph.ac.cell.assemble and aleph.ac.engine.contracts")
    assert rewritten == "aleph.components.incumbent.assemble and aleph.engine.contracts"


def test_a_longer_identifier_is_not_a_match() -> None:
    """`aleph.ac.engineering` would be a different package, and the boundary has to say so."""
    rewritten, count = rewrite("aleph.ac.engineering.thing")
    assert count == 0
    assert rewritten == "aleph.ac.engineering.thing"


def test_the_oracle_leaves_the_law_layer() -> None:
    """`cytosim_parity` is the one check in this tree that is not judged by this tree.

    It goes to `validation/`, where the layer test forbids `laws/` from importing it. An oracle a
    runtime can import is an oracle the runtime can be made to agree with.
    """
    rewritten, _ = rewrite("from aleph.ff.cytosim_parity import compare")
    assert rewritten == "from aleph.validation.cytosim_parity import compare"


def test_the_defective_gamma_floor_path_is_not_archived() -> None:
    """It was, and that was wrong — the modules carry the guards that record their own defects.

    `ff/ENGINE.md`'s banner names three live defects in this path and `STATE.md` (c) 2 blocks every
    number it produces, which is what made archiving look right. But
    `test_three_2026_07_23_defects_are_still_live_in_the_gamma_floor_builder` imports the module,
    so archiving it deletes the record of the defect. **A blocked result and a dead module are
    different things.**
    """
    for module in ("gamma_floor", "gamma_floor_dynamic", "gamma_floor_sweep"):
        rewritten, count = rewrite(f"from aleph.ff.{module} import x")
        assert rewritten == f"from aleph.laws.{module} import x"
        assert count == 1


def test_commons_three_load_bearing_modules_join_the_laws() -> None:
    for module in ("surface_manifold", "turgor_pi0", "compartments"):
        rewritten, _ = rewrite(f"from aleph.common.{module} import x")
        assert rewritten == f"from aleph.laws.{module} import x"


def test_virtual_cell_is_not_moved() -> None:
    """Revised A4: all 26 modules have tests, so it stays until `inner/`/`outer/` can take it."""
    assert not any(old.startswith("virtual_cell") for old, _ in MOVES)
    rewritten, count = rewrite("from aleph.virtual_cell.observation_operator import Op")
    assert count == 0
    assert rewritten == "from aleph.virtual_cell.observation_operator import Op"


def test_nothing_is_rewritten_twice() -> None:
    """Idempotence, or a retry after a partial run corrupts the tree."""
    once, _ = rewrite("from aleph.ac.engine import x\naleph/ff/units.py")
    twice, count = rewrite(once)
    assert twice == once
    assert count == 0


def test_every_move_source_exists_or_has_already_moved() -> None:
    """A typo in MOVES would silently move nothing and rewrite imports to a package with no files."""
    import pathlib

    repo = pathlib.Path(__file__).resolve().parents[3]
    for old, new in MOVES:
        source = repo / "aleph" / pathlib.Path(*old.split("."))
        target = repo / "aleph" / pathlib.Path(*new.split("."))
        assert source.is_dir() or target.is_dir(), f"neither {old} nor {new} exists"
