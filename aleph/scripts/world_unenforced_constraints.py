#!/usr/bin/env python
r"""Constraints a law module STATES that no builder can be given.

**Why this exists.** On 2026-08-21 the microtubule aster was found with 68.6% of its nodes inside the
nuclear envelope, while `aleph/laws/microtubule.py:40` says, as a statement of physical fact, that a
microtubule *"cannot interpenetrate the nucleus"*. The builder has no nucleus argument, so the
sentence could not be honoured even by a caller who had read it. The line directly above it in
`populations.py` calls a builder that DOES take `R_nuc_um`.

**No number in the tree contradicted any other number.** `assert_partitioned` passed — partition is a
property of ID ranges and does not know where a node is. The population counts were right. The render
ran. A defect with no contradiction cannot be caught by computing, which is why it was caught by
looking at a picture (session 21's framing, adopted). This script is the part of that class that CAN
be mechanised: not *"is the constraint satisfied"* — that needs the geometry — but *"could it be"*.

⚠ **THIS IS NOT A GATE AND MAY NEVER BE ONE.** It reports candidates for a person to read. Natural
language is not parsed, only matched: a constraint phrased without one of the declared verbs is
invisible to it, and a builder that honours a constraint through a differently-named argument reads
as a false positive. **Its output is a reading list, and a clean run proves nothing.**

⚠ **And it refuses to be vacuous.** A run that finds no constraint sentences at all, or that fails to
re-find the microtubule/nucleus instance it was written for, EXITS NON-ZERO — because a checker whose
subject came up empty and reported success is the exact failure this whole class is about.

⚠ **MEASURED PRECISION, 2026-08-21: 1 of 15.** The list was read, not just produced. Of fifteen rows
on that day exactly one was actionable — the canary. Eleven are a constraint verb and a subject name
landing in the same sentence of a module docstring that is not a constraint on any builder
(`cortex_assembly:408`, `fa_anchor:1`, `membrane_surface:1`, `turgor_constants:1`, `microtubule:64`,
which describes a defect already fixed). Three are real but already recorded in the module that
raised them — `gamma_floor:154/399/524` on the legacy FF `build_crosslinked_cortex` path, and
`nucleus_envelope:1` saying the FF nucleus *"cannot behave as a coherent"* body because it is a bead
cloud on a radial spring.

⚠ **AND THE RECALL IS NOT MEASURED.** *"Fifteen rows, one real"* reads as *"there is one"*, which is a
recall conclusion drawn from a precision measurement — the same shape as a check whose coverage is
unknown being read as coverage. **This does not know what it missed**, and two limits are structural:
`CONSTRAINT_VERBS` is a fixed list of seven, so *"is prohibited"*, *"shall not"* and *"no X may"* are
invisible; and a constraint must name another builder's subject as a **literal token**, so one written
about *"the organelle"* or *"the surrounding shell"* is invisible too. **A clean run is not evidence
that no unenforceable constraint exists**, and the two warnings the script prints are per-ROW warnings,
not a statement about the set.

**One in fifteen is what a keyword match over prose buys, and it is written here so the next reader
knows the ratio before spending an hour on it.** Improving it would mean parsing the sentences, which
is a different program. **The value is not the precision; it is that the one row it found could not
have been found by comparing any two numbers in the tree.**

Sanity Gate:
    * dimensions — none; this reads source text and function signatures, no physical quantity.
    * boundary cases — zero laws, zero builders, and zero constraint hits are each an explicit refusal
      rather than a pass. A law module with no docstring is reported as unscanned, not as clean.
    * conservation — not applicable.
    * sign sense — not applicable.
    * measurement protocol — every law file gets the same verb list and the same noun extraction, so a
      difference between rows is a difference between modules and not between treatments.

engine units: none. Runtime: host; imports no device module.
"""

from __future__ import annotations

import argparse
import ast
import re
import subprocess
import sys
from pathlib import Path

#: Verbs that make a sentence a CONSTRAINT rather than a description. Declared, not inferred, and
#: deliberately short: a longer list would match prose about what the code does.
CONSTRAINT_VERBS = ("cannot", "can not", "must not", "may not", "never", "is forbidden", "is not allowed")

#: The instance this script was written for. If a run does not re-find it — **naming the nucleus as
#: the subject, not merely containing the word** — the run is a failure. A loose match that happens to
#: catch the canary for the wrong reason is how a checker keeps passing after it stops working.
CANARY = ("microtubule", "nucleus")

#: Subjects the laws constrain that are NOT the stem of any builder module.
#:
#: ⚠ There is exactly one, and it is not an oversight in this script — **it is a gap in the tree, and
#: the tree already knows.** The laws talk about "the nucleus" and the builders have no single name for
#: it: it is assembled from `envelope.py`, `lamina.py` and `chromatin.py`, so a constraint written
#: about the nucleus names nothing any builder is called.
#:
#: **This is NOT a new finding.** `test_families_census.py` already carries it as a strict XFAIL:
#: *"`if_nucleus_linc`'s SPEC names `nucleus`; `build/envelope.py` claims `nuclear_envelope` …
#: Editing either side to make this pass would be choosing which name is right, and that is a PI
#: call."* This line is downstream of that open disagreement, not a discovery of it — and when the PI
#: rules, that XFAIL becomes an XPASS and this list should shrink in the same commit.
#:
#: Without this line the canary matches through the word "filaments" in the same sentence and the row
#: reports the wrong subject — which is how it first ran.
EXTRA_SUBJECTS = ("nucleus",)


def _sentences(text: str) -> list[str]:
    """Split on sentence enders and newlines, keeping comment/docstring prose readable."""
    return [s.strip() for s in re.split(r"(?<=[.;])\s+|\n\s*\n", text) if s.strip()]


def _prose(path: Path) -> list[tuple[int, str]]:
    """Every comment line and docstring in a module, with its 1-indexed line number."""
    src = path.read_text(errors="replace")
    out: list[tuple[int, str]] = []
    for i, line in enumerate(src.splitlines(), 1):
        stripped = line.strip()
        if stripped.startswith("#"):
            out.append((i, stripped.lstrip("# ")))
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return out
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            doc = ast.get_docstring(node)
            if doc:
                out.append((getattr(node, "lineno", 1), doc))
    return out


def main(argv: list[str] | None = None) -> int:
    """Report law-module constraints whose subject is not a parameter of any builder."""
    root = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--laws", type=Path, default=root / "aleph" / "laws")
    ap.add_argument("--builders", type=Path, default=root / "aleph" / "world" / "build")
    args = ap.parse_args(argv)

    law_files = sorted(p for p in args.laws.glob("*.py") if not p.name.startswith("_"))
    build_files = sorted(p for p in args.builders.glob("*.py") if not p.name.startswith("_"))
    if not law_files or not build_files:
        print(f"REFUSED: {len(law_files)} law modules and {len(build_files)} builder modules. "
              "An empty subject is not a clean run.")
        return 2

    # Every parameter name every builder accepts, per module stem. Signatures are read from SOURCE,
    # not by importing -- a builder module imports warp, and this must run where there is no device.
    params: dict[str, set[str]] = {}
    for bf in build_files:
        names: set[str] = set()
        try:
            tree = ast.parse(bf.read_text(errors="replace"))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith("build"):
                a = node.args
                for arg in (*a.posonlyargs, *a.args, *a.kwonlyargs):
                    names.add(arg.arg.lower())
        params[bf.stem] = names

    #: The nouns a constraint can be ABOUT: the names of the things this tree builds, plus the ones
    #: it builds under no single name (see EXTRA_SUBJECTS).
    subjects = {s for s in params if s} | set(EXTRA_SUBJECTS)

    hits: list[tuple[str, int, str, str, bool]] = []
    for lf in law_files:
        stem = lf.stem
        for lineno, block in _prose(lf):
            for sent in _sentences(block):
                low = sent.lower()
                if not any(v in low for v in CONSTRAINT_VERBS):
                    continue
                for other in sorted(subjects):
                    # A constraint in module A that names module B's subject. ⚠ The FULL name must
                    # appear -- matching the last token instead made every mention of "filaments"
                    # a hit on `intermediate_filament`, which is how the canary first passed for the
                    # wrong reason. A narrower match misses more and lies less.
                    noun = other.replace("_", " ")
                    if other == stem or (noun not in low and other not in low):
                        continue
                    key = noun.split()[-1]
                    enforceable = any(key in p or other in p for p in params.get(stem, set()))
                    hits.append((stem, lineno, other, sent, enforceable))

    if not hits:
        print(f"REFUSED: scanned {len(law_files)} law modules and matched NOTHING. Either the verb "
              f"list {CONSTRAINT_VERBS} no longer matches how constraints are written here, or the "
              "prose extraction broke. A checker whose subject came up empty has not passed.")
        return 2

    unenforced = [h for h in hits if not h[4]]
    print(f"scanned {len(law_files)} law modules against {len(build_files)} builders\n")
    print(f"{'law module':<22} {'line':>6}  {'constraint names':<22} enforceable?")
    for stem, lineno, other, sent, ok in hits:
        print(f"{stem:<22} {lineno:>6}  {other:<22} {'yes' if ok else 'NO -- no such parameter'}")
        print(f"        {sent[:150]}")
    print()

    # ⚠ STRICT: the canary must be found as the row's SUBJECT, not merely as a word in the sentence.
    canary_found = any(h[0] == CANARY[0] and h[2] == CANARY[1] and not h[4] for h in hits)
    if not canary_found:
        print(f"REFUSED: the canary is gone. This script exists because {CANARY[0]}'s law module says "
              f"a filament cannot interpenetrate the {CANARY[1]} while no builder takes it as a "
              "parameter. A run that no longer reports that ROW -- subject and all -- is not "
              "reporting on the tree, it is reporting on itself. If the builder was fixed, retire "
              "the canary in the same commit and say so; do not loosen the match.")
        return 2

    print(f"⚠ {len(unenforced)} of {len(hits)} constraints name something no builder in that module "
          "takes as a parameter. That is NOT a verdict -- read them. A constraint honoured through a "
          "differently-named argument reads here as unenforceable and is not.")
    print("⚠ And a clean row proves nothing: this asks whether a constraint COULD be honoured, never "
          "whether it WAS. Whether it was needs the geometry.")
    return 0


def _demo() -> None:
    """Sanity Gate — the parts that can be checked without the tree being any particular way."""
    assert _sentences("a. b\n\nc") == ["a.", "b", "c"]
    # ⚠ PRECONDITION OF THE POSITIVE CONTROL BELOW, not a test of anything. If this ever stops being
    # true, the planted fixture stops exercising the matcher and the control has stopped controlling
    # for anything -- **and this assert fires FIRST, so the control never runs and never reports.**
    # That is exactly what happened on 2026-08-21: breaking CONSTRAINT_VERBS to check the control
    # killed the run here instead. A guard placed before the real test can hide the real test's
    # result; labelling it is the only handling session 21 and this session know of.
    assert any("cannot" in v for v in CONSTRAINT_VERBS), (
        "CONSTRAINT_VERBS no longer contains 'cannot' -- the fixture below says 'a widget cannot "
        "interpenetrate the sprocket', so the positive control has stopped controlling")
    # prose extraction must see both a comment and a docstring
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        f = Path(d) / "m.py"
        f.write_text('"""Doc says it cannot happen."""\n# comment says it must not\nX = 1\n')
        got = " ".join(t for _, t in _prose(f))
        assert "cannot" in got and "must not" in got, got
    # a syntactically broken module must degrade to comments, not raise
    with tempfile.TemporaryDirectory() as d:
        f = Path(d) / "b.py"
        f.write_text("# never do this\ndef (:\n")
        assert _prose(f) == [(1, "never do this")]

    # ⚠ POSITIVE CONTROL, and the reason it is a fixture rather than the tree.
    #
    # `CANARY` is the microtubule/nucleus row -- which is also the ONLY actionable row the scanner has
    # ever produced. So the day PI item 20 is ruled and `build_microtubules` takes a nucleus argument,
    # the canary disappears and this script refuses: **fixing the defect would turn the tool off.**
    # A positive control has to be something that will never be fixed, so it is planted here. This pair
    # is deliberately unenforceable and is not a defect in anything, so it cannot be repaired away.
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "laws").mkdir()
        (Path(d) / "build").mkdir()
        (Path(d) / "laws" / "widget.py").write_text(
            '"""A widget cannot interpenetrate the sprocket."""\n')
        (Path(d) / "build" / "widget.py").write_text("def build_widget(arena, *, seg_um):\n    return {}\n")
        (Path(d) / "build" / "sprocket.py").write_text("def build_sprocket(arena):\n    return {}\n")
        out = subprocess.run(
            [sys.executable, __file__, "--laws", str(Path(d) / "laws"),
             "--builders", str(Path(d) / "build")], capture_output=True, text=True)
        # It refuses on the missing tree canary -- correctly, this is not the tree -- but the ROW must
        # be there, because that row is the whole claim that the matching works.
        assert "widget" in out.stdout and "sprocket" in out.stdout, out.stdout
        assert "NO -- no such parameter" in out.stdout, out.stdout

    print("unenforced-constraints self-check OK — including the planted positive control, which "
          "survives PI item 20 being fixed")


if __name__ == "__main__":
    _demo() if "--self-check" in sys.argv else sys.exit(main())
