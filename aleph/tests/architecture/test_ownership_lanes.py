"""Lane globs must not silently cover the same path, which the loader does not check.

`aleph/coordination/ownership.py` refuses two lanes that declare the **same glob string**, and says
so at load time because *"overlapping ownership is exactly the state this file exists to prevent."*
It does not refuse two *different* globs that resolve to one path — `aleph/engine/**` and
`aleph/engine/cortex_state.py` load without complaint, and then two sessions both believe they own
that file.

That is not hypothetical. It happened the first time a lane was split: on 2026-08-09 a Card-5 session
needed five paths that sat inside `engine`'s globs, and the narrowing had to be expressed as
character-class patterns — `aleph/engine/co[^r]*.py` and friends — which loaded cleanly, read as
noise, and **did not actually exclude anything**. This module is what caught that.

So the rule this file enforces is: **every overlap is either absent or declared here**, and **the
declared winner is the one `owner_of` actually returns.** One overlap is declared, because Card-5's
approval is specific to files `engine` otherwise owns and splitting the whole directory to express it
would read worse than the exception.

How the winner is decided, and why it is worth a test
-----------------------------------------------------
`owner_of` returns on the **first** matching glob and a YAML mapping loads in document order, so the
lane declared first wins. `fnmatch` has no negation — a leading `!` is a literal character — so
exclusion cannot be written at all, and `card5` is placed above `engine` instead.

**Nothing in the file's syntax says the order matters.** Alphabetising the blocks, or adding a lane
at the top, reassigns paths silently. That is the same shape as everything in
`docs/decisions/FINDING-2026-08-09-checks-that-pass-for-the-wrong-reason.md`: position is a proxy for
specificity, and it agrees with specificity until somebody edits the file for an unrelated reason.
"""

from __future__ import annotations

import fnmatch
import itertools
import pathlib

import pytest

yaml = pytest.importorskip("yaml")

ROOT = pathlib.Path(__file__).resolve().parents[3]
OWNERSHIP = ROOT / "ownership.yaml"

#: `(path, winner, loser)` — overlaps that are intended, with the lane that wins by agreement.
#: **Every entry is a place the mechanism cannot express something the project decided**, so each
#: needs a reason in `ownership.yaml` beside the lane, not only here.
#:
#: ⚠ The first two losers changed from `engine` to `archive-readonly` on 2026-08-22, when the PI
#: ruled that the arena is the engine and `aleph/engine/**` moved lanes. This test's own docstring
#: had asked for exactly that — *"if `engine` is narrowed later so it no longer covers these, the
#: exception should go with it"* — and it FAILED on the narrowing rather than passing stale, which
#: is the check working. What the new loser means is sharper than the old one: `card5` may edit two
#: files inside a tree nobody is supposed to edit. That is surfaced beside `archive-readonly` in
#: `ownership.yaml` and left for the PI, because retiring `card5` is retiring a specific approval.
DECLARED_OVERLAPS: tuple[tuple[str, str, str], ...] = (
    ("aleph/engine/cortex_state.py", "card5", "archive-readonly"),
    ("aleph/engine/interior_column_slice.py", "card5", "archive-readonly"),
    ("aleph/scripts/ac_gate_b_interior_column_native.py", "card5", "engine"),
    ("aleph/tests/ac/engine/test_cortex_state.py", "card5", "engine"),
    ("aleph/tests/ac/engine/test_interior_column_ownership.py", "card5", "engine"),
)


def _lanes() -> dict[str, list[str]]:
    return yaml.safe_load(OWNERSHIP.read_text())["sessions"]


def _tracked_paths() -> list[str]:
    import subprocess

    done = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files"], capture_output=True, text=True, check=False
    )
    return done.stdout.split()


def _owners(path: str, lanes: dict[str, list[str]]) -> set[str]:
    """Every lane whose globs cover `path`, matching the way `check_paths` does."""
    hit = set()
    for lane, globs in lanes.items():
        for pattern in globs:
            if fnmatch.fnmatch(path, pattern) or fnmatch.fnmatch(path, pattern.rstrip("/*") + "/*"):
                hit.add(lane)
                break
    return hit


def test_the_declaration_loads() -> None:
    """If it does not load, enforcement is off and every commit passes."""
    from aleph.coordination.ownership import load_ownership

    assert load_ownership(OWNERSHIP) is not None


def test_no_two_lanes_declare_the_same_glob_string() -> None:
    """What the loader already checks — asserted here so a loader regression is visible."""
    lanes = _lanes()
    clashes = [
        (a, b, sorted(set(lanes[a]) & set(lanes[b])))
        for a, b in itertools.combinations(sorted(lanes), 2)
        if set(lanes[a]) & set(lanes[b])
    ]
    assert not clashes, clashes


def test_every_path_level_overlap_is_declared() -> None:
    """The one the loader does NOT check, and the reason this module exists."""
    lanes = _lanes()
    declared = {path for path, _, _ in DECLARED_OVERLAPS}
    undeclared = []
    for path in _tracked_paths():
        if path in declared:
            continue
        owners = _owners(path, lanes)
        if len(owners) > 1:
            undeclared.append(f"{path} → {sorted(owners)}")
    assert not undeclared, (
        "these paths are covered by more than one lane and nobody agreed which wins. Either narrow "
        f"a glob or add the pair to DECLARED_OVERLAPS with its reason: {undeclared[:12]}"
    )


@pytest.mark.parametrize(("path", "winner", "loser"), DECLARED_OVERLAPS)
def test_the_declared_winner_is_the_one_the_loader_actually_returns(
    path: str, winner: str, loser: str
) -> None:
    """The table says who wins; this asks `owner_of`, which is what the hook calls.

    **Resolution is by declaration order.** `owner_of` iterates `sessions.items()` and returns on the
    first matching glob; a YAML mapping loads in document order; so the lane declared first wins.
    Nothing in the YAML syntax says that, which is why moving a block — alphabetising, or adding a
    lane at the top — silently reassigns paths. This test is what turns that into a failure.

    It is also why `ownership.yaml` carries a warning above the lanes and why the real fix is to
    score globs by literal-segment count rather than by position: position is a proxy for
    specificity that holds until somebody edits the file for an unrelated reason.
    """
    from aleph.coordination.ownership import load_ownership

    assert load_ownership(OWNERSHIP).owner_of(path) == winner, (
        f"{path} resolves to something other than {winner!r}. If the lane blocks moved in "
        f"ownership.yaml, the order is load-bearing and that is the cause"
    )


@pytest.mark.parametrize(("path", "winner", "loser"), DECLARED_OVERLAPS)
def test_a_declared_overlap_is_still_really_an_overlap(path: str, winner: str, loser: str) -> None:
    """A declaration that has stopped being true is worse than no declaration.

    If `engine` is narrowed later so it no longer covers these, the exception should go with it —
    otherwise this table becomes a list of things somebody once had to think about.
    """
    owners = _owners(path, _lanes())
    assert winner in owners, f"{winner} no longer claims {path}; the exception is stale"
    assert loser in owners, (
        f"{loser} no longer covers {path}, so this is not an overlap any more — drop the entry"
    )


def test_the_exception_carries_its_reason_where_the_lane_is_declared() -> None:
    """A reason that lives only in a test is a reason the next editor of the YAML will not read."""
    text = OWNERSHIP.read_text()
    assert "card5" in text
    assert "the loader does not notice" in text, (
        "ownership.yaml must say, beside the lane, that the overlap is intended and unenforced"
    )
