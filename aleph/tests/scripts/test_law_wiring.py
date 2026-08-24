"""The wiring survey must not read as encouraging when it is not.

⚠ **WHY THIS EXISTS.** The survey's first run reported **31 WIREABLE** — laws that could be bound with
no acceptance criterion. Twenty-six of those thirty-one carry **no `@wp.kernel` at all**: `units`,
`cell_type`, `architecture_spec`, `ecm_library`, `turgor_constants`, `wlc` are host-side specs,
constants and closed forms. "Could be bound today" was a TRUE sentence about a thing that cannot be
bound in the sense a reader takes it, and it inflated the encouraging number five-fold. The
`NO_KERNEL` verdict and `test_wireable_means_there_is_something_to_launch` are the ratchet on that.
"""

from __future__ import annotations

from aleph.scripts.law_wiring import survey

VERDICTS = {"BOUND", "WIREABLE", "BLOCKED_ON_ACCEPTANCE", "NO_KERNEL"}


def test_every_module_gets_exactly_one_verdict_from_the_closed_set() -> None:
    rows = survey()
    assert len(rows) > 30, len(rows)
    assert {r.verdict for r in rows} <= VERDICTS, {r.verdict for r in rows} - VERDICTS
    assert len({r.module for r in rows}) == len(rows), "a module got two verdicts"


def test_wireable_means_there_is_something_to_launch() -> None:
    """A module with no kernel may never be called WIREABLE. The ratchet on the 31-vs-5 inflation."""
    rows = survey()
    empty = [r.module for r in rows if r.verdict == "WIREABLE" and r.n_kernels == 0]
    assert not empty, (
        f"{empty} were called WIREABLE with no @wp.kernel. There is nothing to bind, and reporting "
        "them as ready reads as progress that does not exist.")


def test_a_bound_module_names_its_importer() -> None:
    """BOUND is an evidenced verdict, not an assumption — the file that binds it is recorded."""
    rows = survey()
    bound = [r for r in rows if r.verdict == "BOUND"]
    assert bound, "nothing is bound, which would itself be the finding"
    for row in bound:
        assert row.importers, row.module
        assert all(i.startswith(("aleph/world", "aleph/scripts/world_")) for i in row.importers), row


def test_a_blocked_module_names_the_state_that_blocks_it() -> None:
    """BLOCKED_ON_ACCEPTANCE must point at the persistent state, so a reader can check the call."""
    rows = survey()
    blocked = [r for r in rows if r.verdict == "BLOCKED_ON_ACCEPTANCE"]
    assert blocked, "nothing blocked — with UndefinedAcceptance shipping, that would be suspicious"
    for row in blocked:
        assert row.state_writes, row.module


def test_the_four_adhesion_modules_are_the_blocked_set() -> None:
    """An anchor, not a total. These four are the adhesion chain and they gate every traction claim.

    ⚠ If this ever fails because the set SHRANK, something was bound — check it was not bound past the
    acceptance criterion. If it GREW, a new kinetic law arrived and its blocker should be stated.
    """
    rows = {r.module: r for r in survey()}
    for stem in ("fa_clutch_warp.py", "fa_ecm.py", "fa_maturation.py", "substrate.py"):
        assert rows[stem].verdict == "BLOCKED_ON_ACCEPTANCE", (stem, rows[stem].verdict)


def test_the_heuristic_is_documented_as_biased_and_under_detecting() -> None:
    """The state test reads kernel PARAMETER NAMES. That limit must travel with the verdict.

    ⚠ It found only four blocked modules while claiming to be biased toward blocking, which is a
    signal it UNDER-detects — `piezo` (a channel that opens) and `nucleus_envelope` (rupture) both
    landed WIREABLE and both carry state by their physics. The caveat is pinned so it cannot quietly
    leave the module while the verdicts stay.
    """
    import inspect

    import aleph.scripts.law_wiring as mod
    source = inspect.getsource(mod)
    assert "heuristic" in source.lower()
    assert "biased toward BLOCKED" in source or "biased toward over-reporting" in source
    assert "WIREABLE` IS NOT `SHOULD BE BOUND" in source or "not mean the law belongs" in source
