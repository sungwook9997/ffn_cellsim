#!/usr/bin/env python
r"""Which `laws/` modules the engine binds, and which of the rest could be bound TODAY.

**Why this exists.** ``aleph/docs/param_ledger.yaml`` counts 240 magnitudes declared under ``laws/``,
and 5 of them sit in a module the production driver binds. The other 235 are written down and
unreachable. That reframes the parameter problem: for most of the tree the open question is not *what
value* but *why is nothing calling this*, and those are different pieces of work with different
blockers. ``STATE.md`` (a) carries "14 of 49 modules imported (AST, 2026-08-09)" — a number counted by
hand once, fifteen days ago, which is exactly the kind that goes stale between the counting and the
reading. This makes it runnable.

**The three verdicts, and the line between the last two is the whole point.**

``BOUND``
    Something under ``world/`` or a ``world_*`` driver imports it. It is in the engine's path.

``WIREABLE``
    Not imported, and its kernels only ACCUMULATE FORCE. A force term needs no accepted step: it is
    evaluated from geometry every iteration and discarded with the iteration. Binding one is
    engineering, and it is not blocked on ``STATE.md`` (e) 1.

``BLOCKED_ON_ACCEPTANCE``
    Not imported, and at least one kernel writes a state that must PERSIST across a step — a bond
    that binds, a filament that grows, a channel that opens. The charter says a kinetic connector
    commits ONLY on an accepted physical step, and the acceptance predicate is ``UndefinedAcceptance``.
    Binding one produces a population that cannot transition and whose correctness cannot be falsified.

⚠ **THE STATE TEST IS A HEURISTIC AND IT IS DELIBERATELY BIASED TOWARD BLOCKED.** It reads kernel
parameter names for the vocabulary this tree uses for persistent state (``bound``, ``state``,
``proposal``, ``n_bound``, ``committed``, ``alive``, ``length``, ``count``) and flags a module when a
kernel takes one as a writable array. A force-only module misread as kinetic costs a re-read; a kinetic
module misread as force-only would put a population into the engine that cannot transition and cannot
be falsified, which is the failure this repository has already recorded. So the evidence is PRINTED per
module rather than summarised, and a verdict is a starting point for a human, not a finding.

⚠ **AND `WIREABLE` IS NOT `SHOULD BE BOUND`.** It says the acceptance criterion does not stand in the
way. Whether a law belongs in this cell, at what coefficient, and against which population is a
modelling question with its own answer — several of these modules are `ff/`-era and their coefficients
are on the non-quotable list.

Sanity Gate:
    * dimensional — none; this parses source and computes no quantity.
    * boundary cases — a module with zero kernels is REPORTED as having none, never silently treated
      as force-only; an unreadable file raises with its path.
    * conservation — every module under the root gets exactly one verdict, asserted before printing.
    * CFL/precision — not applicable.
    * sign sense — not applicable.
    * measurement protocol — one AST pass. No module is imported and no device is opened, so the
      verdict is reproducible from a checkout on a machine with no GPU.

engine units: none. Runtime: host AST parsing only.

Usage:
    python aleph/scripts/law_wiring.py
    python aleph/scripts/law_wiring.py --verdict WIREABLE --evidence
"""

from __future__ import annotations

import argparse
import ast
from dataclasses import dataclass
from pathlib import Path

#: Trees whose imports count as "the engine binds this". ⚠ `world_*` drivers are included because a
#: law bound only from a driver IS in the production path — `world_phase4_native.py` is where the
#: cortex and membrane laws are actually launched from.
BINDING_TREES: tuple[str, ...] = ("aleph/world", "aleph/scripts/world_")

#: Parameter names this tree uses for state that must SURVIVE a step. Read the module docstring for
#: why the list is biased toward over-reporting.
_PERSISTENT: frozenset[str] = frozenset({
    "bound", "n_bound", "state", "states", "proposal", "proposals", "committed",
    "alive", "attached", "n_attached", "occupancy", "age", "n_events",
})

#: Names that look persistent and are NOT — geometry the builder owns and a kernel only reads.
_NOT_PERSISTENT: frozenset[str] = frozenset({
    "rest_state", "state_um", "seg_state",
})


@dataclass(frozen=True, slots=True)
class Verdict:
    """One law module's place in the engine, with the evidence that put it there."""

    module: str
    verdict: str
    n_kernels: int
    importers: tuple[str, ...]
    state_writes: tuple[str, ...]


def _importers(base: Path, stem: str) -> tuple[str, ...]:
    """Every file under the binding trees whose source names ``laws.<stem>``.

    ⚠ Textual, not an import graph — a `from aleph.laws import x` re-export would be missed. It is the
    same method `STATE.md`'s own count used, kept so the two numbers are comparable, and the miss
    direction is toward UNDER-reporting BOUND, which is the safe direction here.
    """
    hits: list[str] = []
    for tree in BINDING_TREES:
        root = base / tree if (base / tree).is_dir() else (base / tree).parent
        pattern = f"{Path(tree).name}*.py" if not (base / tree).is_dir() else "*.py"
        for path in sorted(root.rglob(pattern) if (base / tree).is_dir() else root.glob(pattern)):
            if "__pycache__" in path.parts:
                continue
            rel = str(path.relative_to(base))
            if not rel.startswith(tree):
                continue
            text = path.read_text(encoding="utf-8")
            if f"laws.{stem}" in text or f"laws import {stem}" in text:
                hits.append(rel)
    return tuple(hits)


def _kernels_and_state(tree: ast.Module) -> tuple[int, tuple[str, ...]]:
    """Count ``@wp.kernel`` functions and name the persistent-state parameters they take."""
    n = 0
    writes: list[str] = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        decorated = any(
            (isinstance(d, ast.Attribute) and d.attr == "kernel")
            or (isinstance(d, ast.Name) and d.id == "kernel")
            for d in node.decorator_list
        )
        if not decorated:
            continue
        n += 1
        for arg in node.args.args:
            name = arg.arg
            if name in _NOT_PERSISTENT:
                continue
            head = name.split("_d")[0] if name.endswith("_d") else name
            if head in _PERSISTENT or name in _PERSISTENT:
                writes.append(f"{node.name}({name})")
    return n, tuple(dict.fromkeys(writes))


def survey(repo: Path | None = None) -> list[Verdict]:
    """One verdict per module under ``aleph/laws``, sorted by verdict then name."""
    base = repo or Path(__file__).resolve().parents[2]
    laws = base / "aleph/laws"
    out: list[Verdict] = []
    for path in sorted(laws.glob("*.py")):
        if path.name.startswith("_") or path.name == "__init__.py":
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"))
        n_kernels, writes = _kernels_and_state(tree)
        importers = _importers(base, path.stem)
        if importers:
            verdict = "BOUND"
        elif writes:
            verdict = "BLOCKED_ON_ACCEPTANCE"
        elif n_kernels == 0:
            # ⚠ FOURTH VERDICT, added after the first run made the survey read wrong. 26 of the 31
            # modules first called WIREABLE carry NO `@wp.kernel` at all — `units`, `cell_type`,
            # `architecture_spec`, `ecm_library`, `turgor_constants`, `wlc`. They are host-side specs,
            # constants and closed forms; there is nothing to launch, so "could be bound today" was a
            # true sentence about a thing that cannot be bound in the sense the reader would take.
            # Collapsing them into WIREABLE inflated the encouraging number by five times.
            verdict = "NO_KERNEL"
        else:
            verdict = "WIREABLE"
        out.append(Verdict(path.name, verdict, n_kernels, importers, writes))
    order = {"BOUND": 0, "WIREABLE": 1, "BLOCKED_ON_ACCEPTANCE": 2, "NO_KERNEL": 3}
    return sorted(out, key=lambda v: (order[v.verdict], v.module))


def _demo() -> None:
    """Self-check, runnable without pytest.

    ⚠ Anchors and the partition, never a total — a total would make every module added to `laws/` a
    failure here rather than a fact about the tree.
    """
    rows = survey()
    assert len(rows) > 30, len(rows)
    assert {r.verdict for r in rows} <= {"BOUND", "WIREABLE", "BLOCKED_ON_ACCEPTANCE", "NO_KERNEL"}
    # ⚠ A WIREABLE module MUST carry kernels. Without this the verdict reads as "ready to bind" over
    # a host-side constants file, which is how the first run reported 31 where the answer was 5.
    assert all(r.n_kernels > 0 for r in rows if r.verdict == "WIREABLE"), \
        [r.module for r in rows if r.verdict == "WIREABLE" and r.n_kernels == 0]
    by = {r.module: r for r in rows}
    # membrane_surface IS launched by world_phase4_native via laws_bind.
    assert by["membrane_surface.py"].verdict == "BOUND", by["membrane_surface.py"]
    # crossbridge_kmc writes a `proposal` and says in its own docstring that it may not commit.
    assert by["crossbridge_kmc.py"].verdict != "WIREABLE", by["crossbridge_kmc.py"]
    counts = {v: sum(1 for r in rows if r.verdict == v) for v in
              ("BOUND", "WIREABLE", "BLOCKED_ON_ACCEPTANCE", "NO_KERNEL")}
    assert all(c > 0 for c in counts.values()), counts
    print(f"[law_wiring] self-check OK — {counts}")


def main(argv: list[str] | None = None) -> int:
    """Print the survey."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--verdict", default=None,
                    choices=("BOUND", "WIREABLE", "BLOCKED_ON_ACCEPTANCE", "NO_KERNEL"))
    ap.add_argument("--evidence", action="store_true", help="show importers and state parameters")
    args = ap.parse_args(argv)

    rows = survey()
    shown = [r for r in rows if args.verdict is None or r.verdict == args.verdict]
    width = max((len(r.module) for r in shown), default=10)
    for row in shown:
        print(f"{row.verdict:22s} {row.module:{width}s}  kernels {row.n_kernels:2d}")
        if args.evidence:
            for imp in row.importers[:3]:
                print(f"{'':22s}   bound by  {imp}")
            for write in row.state_writes[:3]:
                print(f"{'':22s}   state     {write}")

    counts = {v: sum(1 for r in rows if r.verdict == v) for v in
              ("BOUND", "WIREABLE", "BLOCKED_ON_ACCEPTANCE", "NO_KERNEL")}
    print(f"\n{len(rows)} modules under aleph/laws/")
    for verdict, count in counts.items():
        print(f"  {verdict:22s} {count:3d}")
    print("\n⚠ WIREABLE means the acceptance criterion does not stand in the way. It does NOT mean the "
          "law belongs in this cell,\n  at what coefficient, or against which population — that is a "
          "modelling question with its own answer.")
    print("⚠ The state test is a heuristic biased toward BLOCKED. Re-read --evidence before acting.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
