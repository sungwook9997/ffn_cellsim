"""`aleph/outer/` imports what exists here, with one quarantine that is named rather than hidden.

53,602 lines arrived on 2026-08-09 from `Project_Aleph`'s `codex/external-training-corpus` branch —
the external training corpus for the outer network. Two kinds of breakage were possible on the way
in and this module holds the line on both.

**The path rename, which was mechanical.** The corpus lived at `corpus/external_training/` and now
lives at `aleph/outer/`, so fourteen files carried `corpus.external_training.…` imports. Rewritten.
`test_no_module_still_imports_the_old_corpus_path` is what stops one coming back with the next
acquisition run.

**The missing dependency, which is a decision.** Six modules under `ragtagcag/` import
`aleph.harness` — `Project_Aleph`'s own knowledge-base substrate, 3,831 lines of CAG, vault, object
store, ledger and snapshot. It was **not** brought, because decision B4 keeps this tree's knowledge
base (Notion contract-graph + DuckDB TAG + Obsidian projection) and `ragtagcag` is a second answer to
the same question — the name is RAG + TAG + CAG.

So those six do not import, deliberately, and `ragtagcag/README.md` says why. The tests below pin the
quarantine at its exact size: **it may not grow, and it may not be quietly resolved by pulling a
second knowledge base in behind it.** If `harness/` is the answer, it is the answer to a port-ledger
entry that argued it, not to an import error.
"""

from __future__ import annotations

import ast
import pathlib

import pytest

OUTER = pathlib.Path(__file__).resolve().parents[2] / "outer"

#: The six `ragtagcag/` modules that import `aleph.harness`. Sorted, relative to `outer/`.
#: **May only shrink**, and only by a port-ledger entry that decides the KB question — not by
#: importing `harness/`.
QUARANTINE: tuple[str, ...] = (
    "ragtagcag/cli.py",
    "ragtagcag/ingest.py",
    "ragtagcag/second_pass/cli.py",
    "ragtagcag/second_pass/combined.py",
    "ragtagcag/second_pass/test_combined.py",
    "ragtagcag/tests/test_ingest.py",
)


def _modules() -> list[pathlib.Path]:
    return [p for p in sorted(OUTER.rglob("*.py")) if "__pycache__" not in p.parts]


def _imports(path: pathlib.Path) -> set[str]:
    try:
        tree = ast.parse(path.read_text(errors="ignore"))
    except SyntaxError:
        return set()
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            found.add(node.module)
    return found


def test_the_corpus_is_actually_here() -> None:
    """A guard over an empty directory is not a guard."""
    modules = _modules()
    assert len(modules) > 200, f"only {len(modules)} modules under outer/; the corpus is 239"
    lines = sum(len(p.read_text(errors="ignore").splitlines()) for p in modules)
    assert lines > 50_000, f"{lines} lines; the corpus is ~53,600"


def test_every_module_parses() -> None:
    """Syntax, before anything else. 645 files were copied across a package boundary."""
    broken = []
    for path in _modules():
        try:
            ast.parse(path.read_text(errors="ignore"))
        except SyntaxError as exc:
            broken.append(f"{path.relative_to(OUTER)}:{exc.lineno}")
    assert not broken, broken


def test_no_module_still_imports_the_old_corpus_path() -> None:
    """`corpus.external_training.*` was the address before the move and resolves to nothing now."""
    offenders = [
        f"{path.relative_to(OUTER)} imports {module}"
        for path in _modules()
        for module in sorted(_imports(path))
        if module.startswith("corpus.")
    ]
    assert not offenders, offenders


def test_the_harness_quarantine_is_exactly_these_six() -> None:
    """It may only shrink, and not by importing a second knowledge base."""
    found = tuple(
        sorted(
            str(path.relative_to(OUTER))
            for path in _modules()
            if any(m == "aleph.harness" or m.startswith("aleph.harness.") for m in _imports(path))
        )
    )
    assert found == QUARANTINE, (
        "the `aleph.harness` quarantine moved. Shrinking it is progress and needs the port-ledger "
        f"entry that decides the KB question; growing it is not. expected {QUARANTINE}, got {found}"
    )


def test_the_quarantine_is_documented_where_somebody_will_hit_it() -> None:
    readme = OUTER / "ragtagcag" / "README.md"
    assert readme.is_file(), "the quarantine needs its reason next to the code, not only in a test"
    assert "B4" in readme.read_text(), "the README must name the decision, not just the symptom"


def test_the_quarantine_is_not_collected() -> None:
    """`testpaths` is `["aleph/tests"]`, so `outer/**/tests/` is outside it — by arrangement.

    If that ever changes, six modules with unresolvable imports start erroring at collection. This
    test fails first and says why, rather than the suite failing later and looking like a defect.
    """
    import tomllib

    pyproject = OUTER.parents[1] / "pyproject.toml"
    paths = tomllib.loads(pyproject.read_text())["tool"]["pytest"]["ini_options"]["testpaths"]
    assert paths == ["aleph/tests"], (
        f"testpaths is {paths}; outer/**/tests/ would now be collected and the quarantined modules "
        "would error at import. Resolve the KB question first, or exclude them explicitly."
    )


#: `viz/`'s four modules that import `aleph.runtime` — the layer decision A2 replaced with
#: `engine/`. Same quarantine shape as `ragtagcag/`, different reason: there is nothing to point
#: them at without carrying a second runtime. Rewriting them against `engine/`'s `CellActor` /
#: `CellWorldTransaction` is a port with a ledger entry. **May only shrink.**
#: **Three, not four.** A `grep -l aleph.runtime` said four and included `publisher.py`, which
#: mentions the module in prose and does not import it. Measured by AST, as everything on this tree
#: eventually has to be — it is the twelfth time in one session that a text scan and the parser gave
#: different structural numbers.
VIZ_QUARANTINE: tuple[str, ...] = (
    "boundary_figure.py",
    "cpu_source.py",
    "scene.py",
)


def test_the_viz_runtime_quarantine_is_exactly_these_three() -> None:
    viz = OUTER.parent / "viz"
    found = tuple(
        sorted(
            p.name
            for p in viz.glob("*.py")
            if any(m == "aleph.runtime" or m.startswith("aleph.runtime.") for m in _imports(p))
        )
    )
    assert found == VIZ_QUARANTINE, (
        "the viz/runtime quarantine moved. Shrinking it means somebody rewrote a module against "
        f"engine/; growing it means a second runtime came in. expected {VIZ_QUARANTINE}, got {found}"
    )


def test_the_sixteen_clean_viz_modules_are_actually_there() -> None:
    """The point of carrying the package was the sixteen, not the four."""
    viz = OUTER.parent / "viz"
    modules = [p for p in viz.glob("*.py") if p.name != "__init__.py"]
    assert len(modules) >= 27, f"{len(modules)} viz modules; 8 from ff/ + 19 from Aleph expected"


@pytest.mark.parametrize(
    "module",
    ["aleph.outer.experiment_factory.schema.validate",
     "aleph.outer.experiment_factory.verification.verify_factory",
     "aleph.infer.fisher",
     "aleph.campaign.rng",
     "aleph.observe.operator",
     "aleph.artifacts.digest"],
)
def test_the_rewritten_paths_actually_import(module: str) -> None:
    """The path rewrite is only right if the result loads. Two of the fourteen, checked live."""
    import importlib

    assert importlib.import_module(module) is not None
