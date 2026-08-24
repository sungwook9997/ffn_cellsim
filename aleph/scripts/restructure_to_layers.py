"""Move the package into its layers: `laws/`, `engine/`, `components/`, `viz/`, `archive/`.

Decision **A2**, as narrowed by §1b and §1c of
`docs/decisions/PROPOSAL-the-merge-and-the-rename.md`. Dry run by default; `--apply` to write.

What this does and does not decide
----------------------------------
The layer order is `units -> state -> laws -> engine -> components -> scripts`, and the one rule the
whole merge rests on is that **`laws/` never imports `engine/` or `components/`**. `ff -> ac` is
already at zero imports and `tests/architecture/test_layer_directions.py` keeps it there; this
script only moves files so the directory names say what the test already enforces.

Three things it deliberately does **not** move, each because a measurement said not to:

* **`ac/cell/` is relocated, not archived** — 243 imports, and `ac/engine` is one of the importers.
  It becomes `components/incumbent/`, which keeps the legal `engine -> components` direction and is
  honest about what it is: the whole-cell assembly the engine still binds while it is strangled.
* **`virtual_cell/` stays where it is** until `inner/` and `outer/` have somewhere to put it. All 26
  of its modules have tests.
* **`ff/` moves wholesale** minus the eight `viz_*`, the external parity oracle, and the three
  `gamma_floor*` modules whose defects `ff/ENGINE.md` names in its own banner.

Both spellings are rewritten — `aleph.ac.engine` and `aleph/ac/engine` — plus the component-wise
quoted form (`ROOT / "aleph" / "ac" / "cell"`) that the *rename* missed and that cost twenty-two
tests. That lesson is why this script handles it from the first run rather than the second.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[2]
PKG = "aleph"

#: `(from, to)` as dotted subpaths under the package. Order matters: longest first, so
#: `ac.cell` is matched before `ac`.
MOVES: tuple[tuple[str, str], ...] = (
    ("ac.engine", "engine"),
    ("ac.cell", "components.incumbent"),
    ("ac.motor", "components.motor"),
    ("ac.ecm", "components.ecm"),
    ("ac.fluid", "components.fluid"),
    ("ac.solid", "components.solid"),
    ("ac.nucleus", "components.nucleus"),
    ("ac.weave", "components.weave"),
    ("ac.emergence", "components.emergence"),
    ("ac.lane_d", "archive.lane_d"),
    ("ff", "laws"),
    # `dcm` was here and is deliberately NOT. It was archived on the strength of "178 imports, all
    # internal", and the layer test caught `laws/relax.py` importing its `implicit_overdamped_step`
    # within a minute of the move: that number says dcm does not look out, not that nothing looks
    # in. Reverting it in the tree and leaving it in this table re-archived it on the next run and
    # split the package in two — the fix has to be in both places or the script undoes it.
)

#: Modules pulled out of `ff/` on the way to `laws/`, because they are not laws.
FF_EXTRACTIONS: tuple[tuple[str, str], ...] = (
    ("ff.cytosim_parity", "validation.cytosim_parity"),
    # `gamma_floor{,_dynamic,_sweep}` were here and are deliberately NOT. Their *numbers* are
    # blocked by `STATE.md` (c) 2, which is what made them look archivable — but the modules carry
    # live guards, and one of them is
    # `test_three_2026_07_23_defects_are_still_live_in_the_gamma_floor_builder`. Archiving the
    # module deletes the test that records the defect, which is the opposite of retiring a claim.
    # A blocked result and a dead module are different things.
    # `architecture_metrics` was here and is deliberately NOT — the same mistake as `dcm` in MOVES,
    # made twice in one session. `laws/ecm_library` calls its `parallel_order_parameter`, and
    # `tests/ac/emergence/test_nematic_oracle` checks that call against an oracle. Reverting the
    # move in the tree without deleting the row here re-archived it on the next run.
)

#: `common/`'s three load-bearing modules join the law layer.
#:
#: `turgor_pi0` reads `Path(__file__).with_name("params_turgor.yaml")`, and moving the module
#: without the YAML broke collection in fourteen test files. Whole-directory moves are safe because
#: `git mv` takes the data with them; **an individually extracted module must be checked for
#: siblings.** `COMMON_SIBLINGS` is that check, written down rather than remembered.
COMMON_KEEPERS: tuple[str, ...] = ("surface_manifold", "turgor_pi0", "compartments")
COMMON_SIBLINGS: tuple[tuple[str, str], ...] = (("common/params_turgor.yaml", "laws/params_turgor.yaml"),)

#: What is left of `common/` once the three keepers leave, by who still imports it.
#:
#: A "common" directory is where modules go when nobody wants to decide which layer owns them, and
#: four of these six are now imported by **nothing at all** — which is what that deferral costs when
#: it runs for a few years. `integrity.py` is the sharpest case: its job is to scan the package for
#: structural defects, and it scans subpackages that were deleted on 2026-07-29.
COMMON_REMAINDER: tuple[tuple[str, str], ...] = (
    ("common.cell_geometry", "laws.cell_geometry"),        # scripts + tests still call it; it is geometry
    ("common.sim_realtime", "dcm.sim_realtime"),           # only `dcm/` imports it, so it belongs there
    ("common.checkpoint", "archive.checkpoint"),           # no importer
    ("common.filament_math", "archive.filament_math"),     # no importer
    ("common.integrity", "archive.integrity"),             # no importer, and scans deleted packages
    ("common.production_policy", "archive.production_policy"),  # no importer
)

SUFFIXES = frozenset({".py", ".md", ".toml", ".cfg", ".yaml", ".yml", ".sh", ".txt", ".ini", ".json"})
EXTENSIONLESS = frozenset({"Makefile", "pre-commit", "pre-push", "Dockerfile"})
SKIP_DIRS = frozenset({".git", "__pycache__", ".pytest_cache", ".ruff_cache", "outputs", "node_modules"})

#: This script names every old path in its own tables, so rewriting it corrupts the mapping.
EXEMPT_FILES = frozenset(
    {
        "scripts/restructure_to_layers.py",
        "tests/scripts/test_restructure_to_layers.py",
        # The rename's controls assert what the *rename* produces (`aleph.ac.engine`). Rewriting
        # them here made them assert this script's output against that script's function, which is
        # two tools quietly breaking each other's tests rather than their code.
        "scripts/rename_package_to_aleph.py",
        "tests/scripts/test_rename_package_to_aleph.py",
    }
)


def _pairs() -> list[tuple[str, str]]:
    """Every rename, longest source first so a prefix never shadows a longer match."""
    pairs = list(FF_EXTRACTIONS) + list(MOVES) + list(COMMON_REMAINDER)
    pairs += [(f"common.{name}", f"laws.{name}") for name in COMMON_KEEPERS]
    return sorted(pairs, key=lambda pair: len(pair[0]), reverse=True)


def rewrite(text: str) -> tuple[str, int]:
    """Rewrite dotted, slashed and component-wise-quoted spellings. Returns text and a count."""
    total = 0
    for old, new in _pairs():
        dotted_old, dotted_new = f"{PKG}.{old}", f"{PKG}.{new}"
        slashed_old = f"{PKG}/" + old.replace(".", "/")
        slashed_new = f"{PKG}/" + new.replace(".", "/")

        # `aleph.ac.engine` -> `aleph.engine`, but not `aleph.ac.engineering`
        text, n = re.subn(rf"{re.escape(dotted_old)}(?![A-Za-z0-9_])", dotted_new, text)
        total += n
        text, n = re.subn(rf"{re.escape(slashed_old)}(?![A-Za-z0-9_])", slashed_new, text)
        total += n

        # `"aleph" / "ac" / "cell"` -> `"aleph" / "components" / "incumbent"`. The rename learned
        # this the hard way: a path built component-wise contains neither a dot nor a slash.
        quoted_old = r'"\s*/\s*"'.join(re.escape(part) for part in ([PKG] + old.split(".")))
        quoted_new = '" / "'.join([PKG] + new.split("."))
        text, n = re.subn(rf'"{quoted_old}"', f'"{quoted_new}"', text)
        total += n
    return text, total


def _git(*args: str, apply: bool) -> None:
    if apply:
        subprocess.run(["git", *args], cwd=REPO, check=True)
    else:
        print(f"  git {' '.join(args)}")


def _ensure_package(directory: pathlib.Path, apply: bool) -> None:
    """A new layer directory needs an `__init__.py` or the move produces an unimportable path."""
    if not apply or directory == REPO / PKG:
        return
    directory.mkdir(parents=True, exist_ok=True)
    marker = directory / "__init__.py"
    if not marker.exists():
        marker.write_text(f'"""Layer package: {directory.name}."""\n')


def move_directories(apply: bool) -> None:
    """`git mv` each subpackage into its layer, then each individually extracted module.

    The extractions are separate because they are *files*, and the first draft of this script moved
    only the directories — so `cytosim_parity`'s imports pointed at `validation/` while the file was
    still sitting in `laws/`. An import rewrite without the matching move is the same defect as a
    move without the rewrite, and it is quieter.
    """
    for old, new in MOVES:
        source = REPO / PKG / pathlib.Path(*old.split("."))
        target = REPO / PKG / pathlib.Path(*new.split("."))
        if not source.exists():
            print(f"  (skip, already moved) {old}")
            continue
        _ensure_package(target.parent, apply)
        _git("mv", str(source.relative_to(REPO)), str(target.relative_to(REPO)), apply=apply)

    file_moves = list(FF_EXTRACTIONS) + list(COMMON_REMAINDER) + [
        (f"common.{name}", f"laws.{name}") for name in COMMON_KEEPERS
    ]
    for old, new in file_moves:
        # `ff/` has already become `laws/` above, so an extraction's source is under its new home.
        # In `--apply` the directory move above has already happened, so an extraction's source is
        # under its new home; in a dry run it has not. Check both, or the dry run reports every
        # extraction as "already moved" — which is a lie a reader would act on.
        parts = old.split(".")
        candidates_ = []
        for first in ({parts[0], "laws"} if parts[0] == "ff" else {parts[0]}):
            candidates_.append(REPO / PKG / pathlib.Path(first, *parts[1:-1]) / f"{parts[-1]}.py")
        source = next((c for c in candidates_ if c.exists()), None)
        target = REPO / PKG / pathlib.Path(*new.split(".")[:-1]) / f"{new.split('.')[-1]}.py"
        if source is None:
            print(f"  (skip, not found) {old}")
            continue
        _ensure_package(target.parent, apply)
        _git("mv", str(source.relative_to(REPO)), str(target.relative_to(REPO)), apply=apply)


def candidates() -> list[pathlib.Path]:
    found: list[pathlib.Path] = []
    for path in REPO.rglob("*"):
        if not path.is_file() or SKIP_DIRS & set(path.parts):
            continue
        if path.suffix not in SUFFIXES and path.name not in EXTENSIONLESS:
            continue
        if any(str(path).endswith(exempt) for exempt in EXEMPT_FILES):
            continue
        found.append(path)
    return sorted(found)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    print("directory moves:")
    move_directories(apply=args.apply)

    touched = total = 0
    for path in candidates():
        try:
            before = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        after, count = rewrite(before)
        if count == 0 or after == before:
            continue
        touched += 1
        total += count
        if args.verbose:
            print(f"  {path.relative_to(REPO)}: {count}")
        if args.apply:
            path.write_text(after, encoding="utf-8")

    for old_rel, new_rel in COMMON_SIBLINGS:
        source, target = REPO / PKG / old_rel, REPO / PKG / new_rel
        if source.exists():
            _ensure_package(target.parent, args.apply)
            _git("mv", str(source.relative_to(REPO)), str(target.relative_to(REPO)), apply=args.apply)

    verb = "rewrote" if args.apply else "would rewrite"
    print(f"\n{verb} {total} references across {touched} files")
    if not args.apply:
        print("dry run — nothing written and nothing moved.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
