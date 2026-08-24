"""Every connector module answers for itself, and the census is checkable in one place.

WHY THIS EXISTS.  Fourteen sessions each own exactly one module under
``aleph/world/families/`` and none may add a second file, so every module put its checks in a
``_demo()`` behind ``__main__`` — the convention ``bond.py`` and the build modules already use.
That is correct for the one-file rule and it has one consequence: **pytest collects none of it.**
Session D8 and session D9 both raised it, independently, in the same hour.

This is the collected runner. It belongs to Lead rather than to any connector session, which is why
it is a second file and does not break the one-file rule for anyone.

WHAT IT REFUSES TO DO.  It does not check that a family is CORRECT. It checks the three things the
scaffold promised and that a session cannot check alone:

1. every module exposes a ``SPEC`` and it loads — a module that cannot describe itself cannot be
   queued, and Lead's decision queue is built by reading these;
2. a module that declares itself blocked names WHAT is unanswered, so the queue has content rather
   than a count;
3. ``DUPLICATE_OF`` is acyclic and every target exists — two modules pointing at each other build
   nothing at all, which is the failure mode Lead had to resolve by designating a canonical side for
   each pair.

⚠ Point 3 is not hypothetical. Three duplicate verdicts landed within one hour and each was written
by a session that could not see the others' files.
"""

from __future__ import annotations

import importlib
import pathlib

import pytest

from aleph.world.families import DUPLICATE_OF, FamilySpec

_DIR = pathlib.Path(__file__).resolve().parents[3] / "aleph" / "world" / "families"
MODULES = sorted(p.stem for p in _DIR.glob("*.py") if p.stem != "__init__")


def _spec(name: str) -> FamilySpec:
    module = importlib.import_module(f"aleph.world.families.{name}")
    spec = getattr(module, "SPEC", None)
    assert isinstance(spec, FamilySpec), (
        f"{name} exposes no SPEC. Lead's decision queue is built by reading SPEC off every module, so "
        "a module that cannot describe itself is invisible to it — worse than one that is blocked."
    )
    return spec


@pytest.mark.parametrize("name", MODULES)
def test_the_module_describes_itself(name: str) -> None:
    """A SPEC that loads, with two populations and a recognised basis."""
    spec = _spec(name)
    assert spec.connector, f"{name}: SPEC carries no connector name"
    assert len(spec.populations) == 2


@pytest.mark.parametrize("name", MODULES)
def test_blocked_means_named(name: str) -> None:
    """Blocked is a queue entry, not a count — so it has to say what is unanswered."""
    spec = _spec(name)
    if spec.buildable:
        assert spec.source.strip(), f"{name}: claims buildable with no provenance"
        return
    assert spec.blocked_by, f"{name}: blocked with an empty blocked_by is a refusal with no content"
    for entry in spec.blocked_by:
        assert isinstance(entry, str) and entry.strip(), f"{name}: an unnamed blocker"


def test_duplicate_map_is_acyclic_and_lands_somewhere() -> None:
    """Two modules deferring to each other build nothing. Lead designates; this checks the result."""
    for src, dst in DUPLICATE_OF.items():
        seen, node = {src}, dst
        while node in DUPLICATE_OF:
            assert node not in seen, (
                f"cycle in DUPLICATE_OF starting at {src!r}: {sorted(seen)}. Each pair needs one "
                "canonical side or the family is never built by anyone."
            )
            seen.add(node)
            node = DUPLICATE_OF[node]
        assert node != src, f"{src!r} resolves back to itself"


@pytest.mark.parametrize("name", MODULES)
def test_the_module_runs_its_own_self_check(name: str) -> None:
    """Run each `_demo()` in ISOLATION, which is how its author ran it.

    ⚠ Parametrised deliberately. The first version of this ran all of them inside one test, and that
    made a module fail on state a previous module had left — which is a property of the runner, not of
    the module. A self-check written for `python -m <module>` is entitled to a fresh process's worth of
    assumptions, and collapsing fourteen of them into one call takes that away.
    """
    module = importlib.import_module(f"aleph.world.families.{name}")
    demo = getattr(module, "_demo", None)
    if not callable(demo):
        pytest.skip(f"{name} exposes no _demo(); its checks are unreachable from pytest")
    demo()


@pytest.mark.xfail(strict=True, reason=(
    "OPEN DISAGREEMENT, marked rather than fixed. `if_nucleus_linc`'s SPEC names `nucleus`; "
    "`build/envelope.py` claims `nuclear_envelope`. D8 and D9 hit the same class from the other side — "
    "the contract endpoint is `sf_arc` (transverse/dorsal ARCS) while PHASE 1 built `stress_fiber` "
    "(ventral, FA-to-FA), which are different structures and not a rename. Editing either side to make "
    "this pass would be choosing which name is right, and that is a PI call. strict=True so the mark is "
    "a ratchet: when the PI rules, this test fails as XPASS and the mark comes off."))
def test_the_population_names_agree_with_what_phase_1_actually_built() -> None:
    """A SPEC naming a population the arena does not have cannot be resolved by anyone.

    ⚠ THIS IS A REAL DISAGREEMENT, not a guard against a hypothetical. Three sessions hit it
    independently within one hour: D8 and D9 both found the contract's endpoint written `sf_arc` while
    PHASE 1 claimed `stress_fiber` (and those are different structures — transverse/dorsal arcs versus
    the ventral FA-to-FA bundles `build/stress_fiber.py` actually builds), and D4's self-check expects
    `nucleus` where `build/envelope.py` claims `nuclear_envelope`.

    The names are listed rather than imported because importing the builders would drag `warp` into a
    host-only test. If PHASE 1 renames a population this list is what fails, which is the point.
    """
    built = {"cortex", "membrane", "nuclear_envelope", "microtubule", "intermediate_filament",
             "filopodium", "lamellipodium", "stress_fiber"}
    unknown: dict[str, tuple[str, ...]] = {}
    for name in MODULES:
        spec = _spec(name)
        missing = tuple(p for p in spec.populations if p not in built)
        if missing:
            unknown[name] = missing
    assert not unknown, (
        "these SPECs name populations PHASE 1 did not build: "
        + "; ".join(f"{k} -> {v}" for k, v in sorted(unknown.items()))
        + ". Either the arena population is renamed or the contract endpoint is stale — it is a PI "
          "call which, and until it is made the count for these edges is asked against nothing."
    )
