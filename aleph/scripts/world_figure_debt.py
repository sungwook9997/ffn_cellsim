#!/usr/bin/env python
r"""Run records with no figure beside them — the charter rule, counted instead of assumed.

The charter says *"Every run writes its own figures beside its record, committed with the data."*
On 2026-08-21 the PHASE 4 tau runs were found to have none, and the obvious next question is whether
that was an exception. **It was not an exception and it is also not the crisis the first count said.**

⚠ **The first version of this count was wrong by 2.5x, and the reason is worth more than the number.**
It looked only in the record's own directory for ``*.png`` and reported **146 of 242**. Two figure
conventions are in use here — ``<dir>/*.png`` and ``<dir>/figs/*.png`` — and it knew about one, so
every directory following the other convention was counted as having nothing. The corrected figure is
**58 of 242**. A checker that knows one of two conventions reports on the convention, not on the tree.

⚠ **A record with no figure is not automatically a defect.** Many of these are retired lanes, raw
sweep members whose figure is the sweep's summary, or diagnostics whose output is a table. This
counts; it does not judge, and it is **not a gate**. Its use is to stop the question being answered
from memory.

Sanity Gate:
    * dimensions — none; this counts files.
    * boundary cases — a directory with records and no figure is listed by name; a tree with NO run
      records at all exits non-zero, because an empty subject is not a clean run.
    * conservation — every record is in exactly one bucket; the two counts must sum to the total, and
      that is asserted rather than trusted.
    * sign sense — not applicable.
    * measurement protocol — a "run record" is a ``.json`` carrying ``device`` or ``kind``; the rule
      is stated here so a change in it is visible as a change in the rule.

engine units: none. Runtime: host.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

#: Where a figure is allowed to live relative to its record. ⚠ BOTH are in use, and knowing only the
#: first is what made this script's own first answer wrong by 2.5x.
FIGURE_GLOBS = ("*.png", "figs/*.png", "figures/*.png")


def figures_for(record: Path) -> list[Path]:
    """Every figure that counts as being 'beside' ``record``, under either convention."""
    return [p for g in FIGURE_GLOBS for p in record.parent.glob(g)]


def main(argv: list[str] | None = None) -> int:
    """Count run records with and without a figure, and name the directories with none."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root", type=Path, nargs="?", default=Path("aleph/outputs/ac"))
    args = ap.parse_args(argv)

    records: list[Path] = []
    for f in sorted(args.root.rglob("*.json")):
        if f.name.endswith(".progress.json"):
            continue
        try:
            d = json.loads(f.read_text())
        except (OSError, ValueError):
            continue
        if isinstance(d, dict) and (d.get("device") or d.get("kind")):
            records.append(f)

    if not records:
        print(f"REFUSED: no run records under {args.root}. An empty subject is not a clean run — "
              "either the path is wrong or the definition of a record above no longer matches.")
        return 2

    without = [f for f in records if not figures_for(f)]
    with_ = len(records) - len(without)
    assert with_ + len(without) == len(records)      # every record in exactly one bucket

    print(f"run records under {args.root}: {len(records)}")
    print(f"  a figure beside them (either convention): {with_:>4}  ({100 * with_ / len(records):.0f}%)")
    print(f"  none:                                     {len(without):>4}  "
          f"({100 * len(without) / len(records):.0f}%)")

    dirs: dict[Path, int] = {}
    for f in without:
        dirs[f.parent] = dirs.get(f.parent, 0) + 1
    if dirs:
        print(f"\ndirectories with records and no figure anywhere ({len(dirs)}):")
        for d, n in sorted(dirs.items(), key=lambda kv: (-kv[1], str(kv[0]))):
            print(f"  {n:>3}  {d}")
    print("\n⚠ NOT a gate and not a verdict. A record with no figure may be a retired lane, a raw "
          "sweep member whose figure is the sweep summary, or a diagnostic whose output is a table. "
          "This exists so the question is not answered from memory.")
    return 0


def _demo() -> None:
    """Sanity Gate — that both conventions are seen, which is the thing that was wrong."""
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        (root / "a").mkdir()
        (root / "b" / "figs").mkdir(parents=True)
        (root / "c").mkdir()
        for sub in "abc":
            (root / sub / "r.json").write_text('{"device": "cuda:0"}')
        (root / "a" / "x.png").write_bytes(b"")            # convention 1
        (root / "b" / "figs" / "x.png").write_bytes(b"")   # convention 2 -- the one that was missed
        assert figures_for(root / "a" / "r.json")
        assert figures_for(root / "b" / "r.json"), "the figs/ convention is invisible again"
        assert not figures_for(root / "c" / "r.json")
        assert main([str(root)]) == 0

    with tempfile.TemporaryDirectory() as d:
        assert main([d]) == 2, "an empty subject must refuse, not report 0 of 0 as clean"
    print("figure-debt self-check OK — both conventions seen, empty subject refused")


if __name__ == "__main__":
    import sys
    if "--self-check" in sys.argv:
        _demo()
    else:
        raise SystemExit(main())
