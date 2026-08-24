"""The import graph has a direction, and this is the test that keeps it.

Why this exists, in one measurement
-----------------------------------
On 2026-08-09 the two candidate trees were compared on exactly this axis:

    ffn_cellsim   ff -> ac              0 imports   (now laws -> engine/components)
    Project_Aleph runtime -> vertical  53 imports   (and vertical -> runtime 33: a cycle)

Aleph's cycle was not designed. It came from four files that each arrived for a good reason —
`law_cases.py` (40 of the 53), a parity harness that must see the components it compares kernels
against; `residency.py`, `resident_world.py` and `resident_cortex.py`, the device-side mirror of a
composition that already existed on the host. Every one of them was a sensible commit. Together
they fused the law library and the assembly layer into one object that cannot be reasoned about in
halves.

**That is the failure mode this test exists for: not a bad decision, but four good ones.** A human
reviewer sees one import in one diff. A counter sees fifty-three.

The direction, and why each edge points the way it does
-------------------------------------------------------
::

    units -> state -> laws -> engine -> components -> scripts

- **A law may not know which component holds it.** The same bending law is used by the cortex, the
  stress fibres and the lamellipodium; the moment it imports one of them, the other two either
  import that one too or grow a copy.
- **The oracles judge the runtime and never the reverse.** A runtime that can import its own oracle
  is a runtime that can be made to agree with it, which is not the same as being right. This one is
  harder here than it was in Aleph, whose `validation/` is a sibling of the package — ours is
  *inside* it, so the rule is enforced by module name rather than by directory nesting.

Renaming
--------
`LAYERS` is the whole configuration. When `ff/` becomes `laws/` and `ac/` splits into `engine/` and
`components/`, edit the table; nothing else here changes.
"""

from __future__ import annotations

import ast
import pathlib

import pytest

PACKAGE = pathlib.Path(__file__).resolve().parents[2]
PACKAGE_NAME = PACKAGE.name

#: `(subpackage, must-not-import)`. Order is the layer order; each entry forbids reaching upward.
LAYERS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("units", ("laws", "engine", "components", "virtual_cell", "scripts", "validation", "archive")),
    ("laws", ("engine", "components", "virtual_cell", "scripts", "validation", "archive")),
    ("engine", ("virtual_cell", "scripts", "validation", "archive")),
    ("components", ("virtual_cell", "scripts", "validation", "archive")),
    # The outer network learns from external research. If it ever imports the engine, it has stopped
    # being independent of it and stage 4 stops being runnable in parallel with stage 1.
    ("outer", ("engine", "components", "laws", "dcm", "scripts", "validation", "archive", "inner")),
    # The inner network learns from THIS engine's own accepted runs; the outer one learns from
    # external research. `CLAUDE.md` states them as two layers with two sources, and the only thing
    # that keeps them two is that neither can reach the other: an inner surrogate that imports
    # `outer` is no longer trained on accepted runs alone, and an outer model that imports `inner`
    # has stopped being independent evidence about the engine and become a function of it.
    #
    # Added 2026-08-20 while `inner/` still holds nothing but its `__init__.py`. That is deliberate
    # and it is the cheap moment: this file's own docstring records that the cycle it was written
    # for came from "four good decisions", each defensible in its own diff. A direction declared
    # before the package opens costs one line; the same direction recovered afterwards cost Aleph
    # a fused runtime.
    #
    # `components` is forbidden and `laws`/`engine` are NOT: a surrogate must be able to name the
    # observables it emulates and the units they carry, and forbidding that would force a copy.
    # What it may not do is reach into a specific part's state.
    ("inner", ("components", "outer", "scripts", "validation", "archive")),
    # `world/` became CANONICAL on 2026-08-20 (`STATE.md` (a)) and `engine/`/`components/` became PORT
    # SOURCE. The direction that makes that real rather than declared: the arena may bind `laws/` — the
    # parity run at `78942fc4` launched two of those kernels over arena-addressed arrays unchanged, which
    # is the whole reason `laws/` was not archived with the rest — but it may NOT reach into the trees it
    # replaces. An arena that imports `engine` or `components` has not replaced them, it has wrapped them,
    # and the 68,968 lines being ported FROM would come along by the back door.
    ("world", ("engine", "components", "dcm", "virtual_cell", "scripts", "validation", "archive")),
)

#: Files a layer OWNS but that do not live inside its package directory, as globs under `aleph/`.
#:
#: Added 2026-08-22, after the PI ruled that the arena is the engine and the rest is archived. The
#: `world` row below had held since 2026-08-20 and `aleph/world/**` was clean — but the guard's scope
#: is the PACKAGE, and two files just outside it were still importing the frozen port source:
#: `tests/world/test_bond.py` took `PROVENANCE` from `engine.forces_manifest` (the very import
#: `b203413d` moved to `units/` and then re-exported, so the old spelling kept working), and
#: `scripts/world_phase1_native.py` took the NVML probe from `components/incumbent`. Neither was
#: hidden and neither was noticed, because the thing that would have noticed was not looking there.
#:
#: A layer's test and its driver are not commentary on the layer, they are how it is exercised and how
#: it is run. An arena whose only native driver imports the tree it replaces has not replaced it.
#:
#: Applied ONLY to the layer-direction test, deliberately. The oracle firewall below keeps the package
#: scope it has always had: a runtime that can see its own oracle can be tuned to agree with it, but a
#: TEST reading an oracle is what oracles are for, and widening both rules at once would forbid that
#: as a side effect of fixing something else.
LAYER_EXTRA_FILES: dict[str, tuple[str, ...]] = {
    "world": ("tests/world/**/*.py", "scripts/world_*.py"),
}


#: Nothing that can run in production may import the oracles. Checked separately from `LAYERS`
#: because it is a different rule with a different reason, and collapsing them would let a future
#: edit to the layer table silently drop the firewall.
RUNTIME_SUBPACKAGES = ("units", "laws", "engine", "components", "coordination", "viz", "world")
#: `world` joined 2026-08-20 when it became CANONICAL. It imports no oracle today and the firewall
#: is cheapest to raise before there is anything to tempt it — the same argument as the `inner` row.

#: `outer/` is above the accept boundary — its inputs are literature and experiment, not simulation —
#: so host arithmetic is allowed there and the oracle firewall is checked separately below rather
#: than by adding it to the runtime list. What it must NOT do is reach into the engine.

#: `dcm/` is a runtime and it **does** import the oracles — two of them, and they are the worst two
#: to import. `test_the_dcm_oracle_breach_does_not_grow` holds the count; the breach itself is
#: reported to the PI rather than fixed here, because deleting those imports changes what `dcm/`
#: computes and that is not a merge decision.
DCM_ORACLE_BREACH = 2


def _imports(path: pathlib.Path) -> set[str]:
    """Absolute module names imported by one file. Relative imports are intra-package by definition."""
    try:
        tree = ast.parse(path.read_text(errors="ignore"))
    except SyntaxError:
        return set()
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            found.add(node.module)
    return found


def _files(subpackage: str) -> list[pathlib.Path]:
    directory = PACKAGE / subpackage
    if not directory.is_dir():
        return []
    return [p for p in sorted(directory.rglob("*.py")) if "__pycache__" not in p.parts]


def _targets(module: str, subpackage: str) -> bool:
    """Does `module` name `subpackage`, spelled absolutely or by the bare-package shorthand?

    Both spellings appear in this tree: `aleph.engine.contracts`, and the bare `ac.engine...`
    that twenty-two scripts used until 2026-07-29. The bare form imports a *second copy* of the
    package, so it has to be recognised here rather than quietly missed.
    """
    # Compare whole dotted components, never string prefixes: `aleph.ac` is a prefix of
    # `aleph.action_log`, and a prefix test would report a package that merely shares three
    # letters. The vacuity control below is the reason this is written out rather than assumed.
    parts = module.split(".")
    if parts[0] == PACKAGE_NAME:
        parts = parts[1:]
    return bool(parts) and parts[0] == subpackage


def _layer_files(layer: str) -> list[pathlib.Path]:
    """The layer's package, PLUS the files it owns outside it (`LAYER_EXTRA_FILES`)."""
    found = set(_files(layer))
    for pattern in LAYER_EXTRA_FILES.get(layer, ()):
        found.update(p for p in PACKAGE.glob(pattern) if "__pycache__" not in p.parts)
    return sorted(found)


def _home(path: pathlib.Path) -> str:
    """The top-level `aleph/` subdirectory a file physically lives in."""
    return path.relative_to(PACKAGE).parts[0]


@pytest.mark.parametrize(("layer", "forbidden"), LAYERS, ids=[layer for layer, _ in LAYERS])
def test_a_layer_never_imports_upward(layer: str, forbidden: tuple[str, ...]) -> None:
    violations: list[str] = []
    for path in _layer_files(layer):
        # A file is never forbidden from the directory it lives in. This is a no-op for the package
        # itself — `world/` is not in `world/`'s own forbidden set — and it is what makes an extra
        # file's inherited rule mean the right thing: "the arena package must not depend on drivers"
        # is a statement about `aleph/world/**`, and reading it as "a driver must not import a
        # driver" would forbid `world_phase4_native.py` from calling `run_provenance`, which is the
        # module that stamps a run with the build it ran. Encoded as one condition rather than a
        # carve-out list, so a future extra glob cannot need a new exception.
        home = _home(path)
        for module in sorted(_imports(path)):
            for target in forbidden:
                if target != home and _targets(module, target):
                    violations.append(f"{path.relative_to(PACKAGE.parent)} imports {module}")
    assert not violations, (
        f"`{layer}/` reached upward. Aleph's runtime/vertical cycle reached 53 edges one good "
        f"reason at a time; this is what stops the first one: {violations}"
    )


@pytest.mark.parametrize("subpackage", RUNTIME_SUBPACKAGES)
def test_no_runtime_path_imports_the_oracles(subpackage: str) -> None:
    """A runtime that can see its own oracle can be made to agree with it."""
    violations = [
        f"{path.relative_to(PACKAGE.parent)} imports {module}"
        for path in _files(subpackage)
        for module in sorted(_imports(path))
        if _targets(module, "validation")
    ]
    assert not violations, violations


# ── vacuity controls ──────────────────────────────────────────────────────────────────────────
def test_every_runtime_subpackage_actually_exists() -> None:
    """`common` sat in this list after the directory was deleted, scanning nothing and passing.

    A firewall test whose subject does not exist is indistinguishable from one that holds.
    """
    missing = [name for name in RUNTIME_SUBPACKAGES if not _files(name)]
    assert not missing, f"{missing} named in RUNTIME_SUBPACKAGES but hold no Python files"


def test_every_extra_layer_glob_matches_something() -> None:
    """A widened scope that resolves to no files is indistinguishable from one that holds.

    This is the `common` lesson one row down, applied before it can happen: that name sat in
    `RUNTIME_SUBPACKAGES` after its directory was deleted, scanning nothing and passing green.
    A glob is easier to break than a directory name — a rename of one driver is enough.
    """
    empty = [
        f"{layer}: {pattern}"
        for layer, patterns in LAYER_EXTRA_FILES.items()
        for pattern in patterns
        if not [p for p in PACKAGE.glob(pattern) if "__pycache__" not in p.parts]
    ]
    assert not empty, f"LAYER_EXTRA_FILES globs match no files, so they guard nothing: {empty}"


def test_the_dcm_oracle_breach_does_not_grow() -> None:
    """`dcm/` imports two closed-form oracles, and the charter forbids exactly this.

    *"A closed form (Bell-Evans, Hill, Pereverzev) is an acceptance ORACLE, never the runtime
    mechanism. The v1 codebase used published models AS the runtime; this one inverts that."*

    `dcm/dcm_ecm_clutch_host.py` imports `validation.pereverzev`, and `dcm/dcm_cadherin_host.py`
    imports `validation.cadherin_sliding_rebinding`. Both are the runtime consuming its own oracle —
    the v1 pattern the project was rebuilt to invert, still live in the parked tree.

    Ratcheted rather than fixed: removing the imports changes what `dcm/` computes, which is a
    physics decision and not a merge one. It may only shrink.
    """
    breaches = [
        f"{path.relative_to(PACKAGE.parent)} imports {module}"
        for path in _files("dcm")
        for module in sorted(_imports(path))
        if _targets(module, "validation")
    ]
    assert len(breaches) <= DCM_ORACLE_BREACH, (
        f"the runtime-imports-its-own-oracle breach grew from {DCM_ORACLE_BREACH} to "
        f"{len(breaches)}: {breaches}"
    )


def test_the_scanner_finds_the_imports_that_are_there() -> None:
    """`ac/` really does import `ff/`. If this drops to zero the scanner has stopped working.

    The permitted direction is the control for the forbidden one: a scanner that returned nothing
    would pass every test above while checking nothing at all.
    """
    downward = [
        path
        for path in _files("components") + _files("engine")
        if any(_targets(module, "laws") for module in _imports(path))
    ]
    assert len(downward) > 20, (
        f"only {len(downward)} files above laws/ import it — measured at 55 imports on 2026-08-09 "
        "(components->laws 40, engine->laws 15), so either the scanner is broken or the engine has "
        "been restructured without this test"
    )


def test_the_target_matcher_distinguishes_a_prefix_from_a_package() -> None:
    assert _targets("aleph.engine.contracts", "engine")
    assert _targets("engine.contracts", "engine")
    assert _targets("aleph.engine", "engine")
    # A different package that merely starts with the same letters is not a violation.
    assert not _targets("aleph.engineering_log", "engine")
    assert not _targets("engines", "engine")
    assert not _targets("aleph.laws.units", "engine")


def test_every_named_layer_actually_exists() -> None:
    """A typo in `LAYERS` would make its row scan an empty directory and pass silently."""
    missing = [layer for layer, _ in LAYERS if not _files(layer)]
    assert not missing, f"{missing} named in LAYERS but hold no Python files"
