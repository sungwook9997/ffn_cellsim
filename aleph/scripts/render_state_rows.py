#!/usr/bin/env python
r"""Render STATE.md's tier-(a) table from an append-only YAML, so parallel sessions stop colliding.

THE PROBLEM THIS SOLVES IS NOT MERGE CONFLICTS. Git reports those. It is the silent case: two
sessions editing the rendered markdown table minutes apart, where the loser's row is not conflicted
but simply absent, or is swept into the winner's commit and attributed to it. Both happened on
2026-07-28 within an hour — one of them carried away an entire commit's staged content.

A markdown table is the worst possible shared format for concurrent writers: every row is on one
physical line, every session appends at a different place, and the surrounding prose shifts. A YAML
LIST is the best available one: appends land at the end, each entry is several lines so a conflict is
localised and legible, and the fields are named so a half-merged entry is visibly broken rather than
quietly wrong.

    python aleph/scripts/render_state_rows.py            # rewrite STATE.md's table in place
    python aleph/scripts/render_state_rows.py --check    # CI: fail if the table is stale
    python aleph/scripts/render_state_rows.py --extract  # one-off: seed the YAML from STATE.md

WHAT IS NOT AUTOMATED, deliberately. The rows are not derived from run records. A tier-(a) row is a
JUDGMENT — "only what the left column says is quotable" — and deriving it from an artifact would
turn every green run into a quotable claim, which is the exact failure the 2026-07-28 retraction was.
The YAML holds the judgment; this script only renders it, and enforces the two structural rules that
can be checked mechanically: every row names a build commit, and every row names its population.

Sanity Gate:
    * boundary: a missing or empty YAML leaves STATE.md untouched and says so, rather than rendering
      an empty table over a populated one.
    * conservation/invariant: `--check` compares rendered-vs-on-disk exactly, so a hand-edit of the
      generated block is caught instead of being silently overwritten later.
    * measurement-protocol: `population` is REQUIRED on every row. The 2026-07-28 retraction turned on
      nobody seeing at a glance that four headline numbers came from 0.18% of native.
    * numerical / dimensional / sign-sense: not applicable — no quantity is computed here.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

BEGIN = "<!-- STATE-TIER-A:BEGIN -->"
END = "<!-- STATE-TIER-A:END -->"

#: Append-only source of truth for the rendered table. In `free:` in the ownership declaration
#: precisely because every session must be able to append to it.
DEFAULT_ROWS = Path("aleph/docs/state_rows.yaml")

HEADER = (
    "| Result — only what the left column says is quotable | Population / device | Build commit "
    "| Committed artifact |\n|---|---|---|---|"
)


def render(rows: list[dict]) -> str:
    """Return the markdown table body for ``rows``.

    Args:
        rows: Entries with ``headline``, ``claim``, ``population``, ``build``, ``artifact``.
            Only ``headline`` is rendered into STATE.md; the full ``claim`` stays in the YAML. The
            split exists because the claim column was 62% of STATE.md's tier-(a) table and 30% of the
            whole file, and it is evidence — read when a row is challenged, not on every boot.

    Returns:
        The table, header included, without the surrounding markers.

    Raises:
        ValueError: If a row omits a required field. Both structural rules live here: a row without a
            build commit cannot be traced to source, and one without a population cannot be read for
            scale — and scale is what the 07-25 audit and the 07-28 retraction each turned on.
    """
    lines = [HEADER]
    for index, row in enumerate(rows):
        missing = [k for k in ("headline", "claim", "population", "build", "artifact")
                   if not str(row.get(k, "")).strip()]
        if missing:
            raise ValueError(
                f"tier-(a) row {index} is missing {missing}. Every row names a build commit AND its "
                f"population — a claim that states neither is what (c) exists to hold instead."
            )
        lines.append(
            f"| {row['headline'].strip()} | {row['population'].strip()} | {row['build'].strip()} "
            f"| {row['artifact'].strip()} |"
        )
    return "\n".join(lines)


def splice(state_text: str, table: str) -> str:
    """Replace the marked block in ``state_text`` with ``table``.

    Args:
        state_text: Current STATE.md contents.
        table: Rendered table.

    Returns:
        The new file contents.

    Raises:
        ValueError: If the markers are absent or out of order — better to refuse than to guess where
            a generated block belongs in a hand-written file.
    """
    if BEGIN not in state_text or END not in state_text:
        raise ValueError(f"STATE.md is missing {BEGIN} / {END}")
    head, rest = state_text.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    return f"{head}{BEGIN}\n\n{table}\n\n{END}{tail}"


def extract(state_text: str) -> list[dict]:
    """Recover rows from an already-rendered STATE.md, to seed the YAML once.

    Args:
        state_text: Current STATE.md contents.

    Returns:
        One dict per data row.
    """
    body = state_text.split(BEGIN, 1)[1].split(END, 1)[0]
    rows = []
    for line in body.strip().splitlines():
        line = line.strip()
        if not line.startswith("|") or line.startswith("|---") or "only what the left column" in line:
            continue
        cells = [c.strip() for c in line.strip("|").split(" | ")]
        if len(cells) == 4:
            rows.append(dict(zip(("headline", "population", "build", "artifact"), cells)))
    return rows


def main() -> int:
    """Render, check, or extract."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--state", default="STATE.md")
    ap.add_argument("--rows", default=str(DEFAULT_ROWS))
    ap.add_argument("--check", action="store_true", help="fail if STATE.md is stale; write nothing")
    ap.add_argument("--extract", action="store_true", help="seed the YAML from the rendered table")
    args = ap.parse_args()

    state_path, rows_path = Path(args.state), Path(args.rows)
    state_text = state_path.read_text()

    if args.extract:
        rows = extract(state_text)
        rows_path.parent.mkdir(parents=True, exist_ok=True)
        rows_path.write_text(
            "# Append-only source for STATE.md's tier-(a) table. Render with `make state`.\n"
            "#\n"
            "# APPEND HERE, never to the rendered table: a markdown table is one physical line per\n"
            "# row and every concurrent session appends at a different place, so the loser's row goes\n"
            "# missing without conflicting. A YAML list conflicts legibly or not at all.\n"
            "#\n"
            "# EVERY ROW NAMES A BUILD COMMIT AND ITS POPULATION. Both are enforced by the renderer.\n"
            "# A claim that can state neither belongs in STATE.md (c), not here.\n\n"
            + yaml.safe_dump({"rows": rows}, allow_unicode=True, sort_keys=False, width=10_000)
        )
        print(f"seeded {rows_path} with {len(rows)} row(s)")
        return 0

    if not rows_path.exists():
        print(f"{rows_path} does not exist — STATE.md left untouched. Seed it with --extract.",
              file=sys.stderr)
        return 1
    rows = (yaml.safe_load(rows_path.read_text()) or {}).get("rows") or []
    if not rows:
        print(f"{rows_path} declares no rows — refusing to render an empty table over a populated "
              f"one", file=sys.stderr)
        return 1

    rendered = splice(state_text, render(rows))
    if args.check:
        if rendered == state_text:
            print(f"[state-rows] OK — {len(rows)} row(s), STATE.md matches {rows_path}")
            return 0
        print(f"[state-rows] STALE — STATE.md's tier-(a) block does not match {rows_path}. "
              f"Run `make state`. (If you hand-edited the block, move the edit into the YAML: the "
              f"block is generated and the next render would silently discard it.)", file=sys.stderr)
        return 1

    state_path.write_text(rendered)
    print(f"[state-rows] rendered {len(rows)} row(s) into {state_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
