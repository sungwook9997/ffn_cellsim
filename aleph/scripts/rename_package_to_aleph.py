"""Rename the package `ffn_sim` -> `aleph`, and refuse to touch the three things that are not it.

Decision **A1** of `docs/decisions/PROPOSAL-the-merge-and-the-rename.md`. Run with `--dry-run`
(the default) until the output is what you expect, then `--apply`.

The hazard this script exists for, measured before it was written
-----------------------------------------------------------------
`ffn_sim` names four different things in this tree, and a plain
``sed -i 's/ffn_sim/aleph/g'`` corrupts three of them:

======================================  =========  ==================================
what                                      files    rename?
======================================  =========  ==================================
the Python package / import path            693    **yes** — this is the rename
the **conda environment**                     38    **no** — `~/miniconda3/envs/ffn_sim`
                                                    is a directory on disk that this
                                                    change does not move
the GitHub remote `sungwook9997/ffn_cellsim`   1    **no** — a remote is renamed on
                                                    the server or not at all
the string inside committed run artifacts   many    **no** — an artifact records what
                                                    ran, and editing it makes the
                                                    record disagree with history
======================================  =========  ==================================

38 files reference the interpreter as `~/miniconda3/envs/ffn_sim/bin/python`. Rewriting those gives
a path that does not exist, and the failure surfaces only when a driver is next launched — on the
GPU host, at the start of a run that then does not happen.

So the substitution is **not** on the bare token. It is on the token in the two positions that mean
"the package": `ffn_sim.` as a dotted module prefix, and `ffn_sim/` as a path segment — with the
environment spellings masked out first and restored afterwards.

What this script does not do
----------------------------
It does not `git mv` the directory, does not commit, and does not touch `ffn_sim/outputs/`. The
directory move is one command a human should run and see; the artifacts are history.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parents[2]
OLD, NEW = "ffn_sim", "aleph"

#: Directories never rewritten. `outputs/` holds committed run records: an artifact that says it was
#: produced by `ffn_sim.scripts.x` is a true statement about a run that happened, and rewriting it
#: would make the record disagree with the commit it names.
SKIP_DIRS = frozenset({".git", "__pycache__", ".pytest_cache", ".ruff_cache", "outputs", "node_modules"})

#: Paths under a skipped directory that are **not** run records and must be rewritten anyway.
#:
#: `outputs/tag_kb/` is the knowledge-base subsystem — its scripts, and the `params_manifest` /
#: `results_manifest` the pre-commit gates read. Skipping it left every manifest pointing at
#: `ffn_sim/configs/...` and `ffn_sim/outputs/...`, so both KB drift gates failed with 66 entries
#: reported as "constant ABSENT on disk". The files had not gone anywhere; the manifest had.
#:
#: That it lives under `outputs/` at all is the reason the structure plan lifts it to a first-class
#: `kb/` package — a live gate configuration filed among immutable run records will be mistaken for
#: one again.
SKIP_EXCEPTIONS: tuple[str, ...] = ("outputs/tag_kb/",)

SUFFIXES = frozenset({".py", ".md", ".toml", ".cfg", ".yaml", ".yml", ".sh", ".txt", ".ini", ".json"})

#: Extensionless files that are nonetheless full of package paths. A suffix allowlist misses these
#: silently, and the first run of this script did: `Makefile` kept nine stale paths and
#: `.githooks/pre-commit` five, so `make kb-check` and every commit would have failed on a path that
#: no longer exists. Named individually rather than by "has no suffix", which would sweep binaries.
EXTENSIONLESS: frozenset[str] = frozenset({"Makefile", "pre-commit", "pre-push", "Dockerfile"})

#: This script and its controls talk *about* the old name, so a rewrite corrupts them.
#:
#: Found by running it: the docstring above became ``sed -i 's/aleph/aleph/g'``, and every case in
#: the control table became ``("from aleph.ac import x", "from aleph.ac import x")`` — a test that
#: asserts a no-op and passes. **A self-referential rewrite silently disarms its own tests**, which
#: is worse than corrupting them, because the suite stays green.
EXEMPT_FILES: frozenset[str] = frozenset(
    {
        "scripts/rename_package_to_aleph.py",
        "tests/scripts/test_rename_package_to_aleph.py",
    }
)

#: Spellings where `ffn_sim` is the **conda environment**, not the package. Masked before
#: substitution and restored after, so the two can never be confused by a regex that grew.
ENVIRONMENT_FORMS: tuple[str, ...] = (
    "envs/ffn_sim",
    "activate ffn_sim",
    "conda activate ffn_sim",
    "ffn_sim/bin/python",
    "-n ffn_sim",
    "name: ffn_sim",
)

#: The remote. Renaming a repository is a server-side action; this is only a string here.
REMOTE_FORMS: tuple[str, ...] = ("sungwook9997/ffn_cellsim", "ffn_cellsim.git")

_MASK = "\x00MASK{}\x00"

#: `ffn_sim` in a position that means *the package*:
#:
#: * followed by `.` or `/` — a dotted module prefix or a path segment;
#: * preceded by `import ` or `from ` — the bare-module spellings, `import ffn_sim` and
#:   `from ffn_sim import x`, which carry no dot and no slash and must still move.
#:
#: A bare mention with neither is prose about the old name — *"the tree was called ffn_sim until
#: 2026-08-09"* — and rewriting it would edit the history the rename is supposed to leave readable.
#: * a **quoted standalone token** — `ROOT / "ffn_sim" / "docs"`. This is the position the first run
#:   missed and it broke twenty-two tests: a path assembled component-wise never contains a dot or a
#:   slash, so nothing about the string says "package" except that a human knows it is one.
#:   `conftest.py`'s `_CANONICAL_ROOT = "ffn_sim"` was the first symptom and the import failed
#:   immediately; the twenty-two that did not fail immediately were quieter.
_PACKAGE = re.compile(
    rf"(?<=\bimport ){OLD}\b"
    rf"|(?<=\bfrom ){OLD}\b"
    rf"|(?<=[\"']){OLD}(?=[\"'*])"
    rf"|\b{OLD}(?=[./])"
)


def _mask(text: str) -> tuple[str, list[str]]:
    """Replace every protected spelling with a sentinel, longest first so prefixes cannot shadow."""
    saved: list[str] = []
    for form in sorted(ENVIRONMENT_FORMS + REMOTE_FORMS, key=len, reverse=True):
        while form in text:
            text = text.replace(form, _MASK.format(len(saved)), 1)
            saved.append(form)
    return text, saved


def _unmask(text: str, saved: list[str]) -> str:
    for index, form in enumerate(saved):
        text = text.replace(_MASK.format(index), form)
    return text


def rewrite(text: str) -> tuple[str, int]:
    """Return the rewritten text and how many package references changed."""
    masked, saved = _mask(text)
    rewritten, count = _PACKAGE.subn(NEW, masked)
    return _unmask(rewritten, saved), count


def candidates() -> list[pathlib.Path]:
    found: list[pathlib.Path] = []
    for path in REPO.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix not in SUFFIXES and path.name not in EXTENSIONLESS:
            continue
        relative = str(path.relative_to(REPO))
        if SKIP_DIRS & set(path.parts) and not any(x in relative for x in SKIP_EXCEPTIONS):
            continue
        if any(str(path).endswith(exempt) for exempt in EXEMPT_FILES):
            continue
        found.append(path)
    return sorted(found)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true", help="write the changes (default is a dry run)")
    parser.add_argument("--verbose", action="store_true", help="list every file touched")
    args = parser.parse_args(argv)

    touched = total = protected = 0
    for path in candidates():
        try:
            before = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        if OLD not in before:
            continue
        after, count = rewrite(before)
        if any(form in before for form in ENVIRONMENT_FORMS):
            protected += 1
        if count == 0 or after == before:
            continue
        touched += 1
        total += count
        if args.verbose:
            print(f"  {path.relative_to(REPO)}: {count}")
        if args.apply:
            path.write_text(after, encoding="utf-8")

    verb = "rewrote" if args.apply else "would rewrite"
    print(f"{verb} {total} package references across {touched} files")
    print(f"protected the conda environment spelling in {protected} files")
    if not args.apply:
        print("\ndry run — nothing written. Re-run with --apply, then:")
        print(f"  git mv {OLD} {NEW}")
        print("  python -m pytest -q")
    return 0


if __name__ == "__main__":
    sys.exit(main())
