#!/usr/bin/env python
r"""Refuse a commit that touches paths this session has not claimed.

Wire it as the pre-commit hook (`make hooks`) or run it by hand:

    FFN_SESSION=sf-motor python aleph/scripts/check_ownership.py
    FFN_SESSION=sf-motor python aleph/scripts/check_ownership.py --paths a.py b.py

With `FFN_SESSION` unset the check is OFF and says so. That is deliberate: this is a guard for
deliberate parallel work, not a tax on one session working alone, and a guard that fires when nobody
asked for it gets disabled and then does not fire when it matters.

What it is actually for is NOT merge conflicts — git reports those. It is the silent case: two
sessions editing one shared status file minutes apart, where the loser's edit is not conflicted but
simply gone, or is swept into the winner's commit. Both happened on 2026-07-28 within an hour, once
carrying away a whole commit's staged content.

Sanity Gate:
    * boundary: no staged paths -> pass; unset session -> pass with a notice, never a silent pass.
    * conservation/invariant: every offending path is reported in one run, so the rule is learned
      once rather than rediscovered file by file.
    * measurement-protocol: paths come from `git diff --cached --name-only`, i.e. exactly what the
      commit will contain — not the working tree, which is a different set.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

# Runs as a git hook, where there is no path setup at all. Import via `aleph.…`, never the bare
# `ac.…` identity that double-loads modules outside pytest.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from aleph.coordination.ownership import (  # noqa: E402
    DEFAULT_OWNERSHIP_PATH,
    OwnershipViolation,
    check_paths,
    load_ownership,
)

SESSION_ENV = "FFN_SESSION"


def _staged_paths(repo: Path) -> list[str]:
    """Return repo-relative POSIX paths staged for commit."""
    done = subprocess.run(["git", "-C", str(repo), "diff", "--cached", "--name-only"],
                          capture_output=True, text=True, check=False)
    return [line.strip() for line in done.stdout.splitlines() if line.strip()]


def main() -> int:
    """Check the staged set against the claim declaration."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--paths", nargs="*", default=None,
                    help="paths to check instead of the staged set")
    ap.add_argument("--ownership", default=None, help="claim file (default ownership.yaml)")
    ap.add_argument("--session", default=None, help="override $FFN_SESSION")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[2]
    session = args.session or os.environ.get(SESSION_ENV, "")
    if not session:
        print(f"[ownership] {SESSION_ENV} unset — parallel-session enforcement is OFF for this commit")
        return 0

    path = Path(args.ownership) if args.ownership else repo / DEFAULT_OWNERSHIP_PATH
    try:
        ownership = load_ownership(path)
    except FileNotFoundError:
        print(f"[ownership] no claim file at {path} — enforcement OFF")
        return 0
    except OwnershipViolation as bad:
        print(f"[ownership] the claim file itself is invalid: {bad}", file=sys.stderr)
        return 1

    paths = args.paths if args.paths is not None else _staged_paths(repo)
    if not paths:
        print("[ownership] nothing staged")
        return 0

    complaints = check_paths(session, paths, ownership)
    if not complaints:
        print(f"[ownership] OK — {len(paths)} path(s), all owned by {session!r}")
        return 0

    print(f"[ownership] REFUSED — session {session!r} may not commit {len(complaints)} of "
          f"{len(paths)} staged path(s):", file=sys.stderr)
    for c in complaints:
        print(f"  - {c}", file=sys.stderr)
    print("\n  Unset FFN_SESSION to bypass (and then you own the merge).", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
