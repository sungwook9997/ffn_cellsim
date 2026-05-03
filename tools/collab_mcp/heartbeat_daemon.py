"""Heartbeat daemon for the ACS Collab Workroom sidebar.

Every 5 s, refresh the `agent_status` row for each of the four LLM panes
in the workroom so the sidebar "Agents" panel stays visibly alive even
when the pane is idle between PI prompts. The daemon never invents
progress numbers — it re-POSTs the existing row's percent/activity
unchanged, which only bumps ``updated_at`` and resets the freshness tick.
If a prefixed workroom row does not exist yet, the daemon seeds a neutral
``0% idle`` row so newly cold-started rooms show their own cards instead
of falling back to stale legacy logical rows.

When ``--workroom <name>`` is given, agent IDs and tmux targets are both
prefixed with the workroom name (e.g. ``design-discussion-claude-chat``).
This lets multiple workrooms heartbeat in parallel without colliding on
the same ``agent_status`` row.

If the underlying tmux pane is dead, the daemon overwrites the row's
activity with ``"(pane dead)"`` and leaves the percent intact so PI can
spot a crashed CLI from the sidebar without inventing fake progress.
"""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import re
import sqlite3
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from contextlib import closing

DB_PATH = pathlib.Path(
    os.environ.get(
        "COLLAB_MCP_DB",
        str(pathlib.Path.home() / ".acs-collab" / "inbox.db"),
    )
).expanduser()
ROOM_URL = os.environ.get("COLLAB_ROOM_URL", "http://127.0.0.1:7879/agent_status")
ROOM_TOKEN = os.environ.get("COLLAB_ROOM_TOKEN", "acs-room")
INTERVAL_S = float(os.environ.get("COLLAB_HEARTBEAT_INTERVAL", "5"))
WORKROOM_RE = re.compile(r"^[a-z0-9][a-z0-9-]{0,31}$")

LOGICAL_PANES = ("claude-chat", "claude-work", "codex-chat", "codex-work")


def build_panes(workroom: str | None) -> list[tuple[str, str]]:
    """Return [(agent_id, tmux_target), ...] for the daemon to refresh.

    Without a workroom (legacy mode), agent IDs and tmux targets are the
    bare logical names (``claude-chat`` etc.). With a workroom, both are
    prefixed (``design-discussion-claude-chat``) so per-room heartbeats
    do not collide on the shared ``agent_status`` table.
    """
    if not workroom:
        return [(name, f"{name}:0.0") for name in LOGICAL_PANES]
    if not WORKROOM_RE.match(workroom):
        raise ValueError(
            f"invalid workroom name {workroom!r}: must match {WORKROOM_RE.pattern}"
        )
    return [(f"{workroom}-{name}", f"{workroom}-{name}:0.0") for name in LOGICAL_PANES]


def pane_alive(target: str) -> bool:
    res = subprocess.run(
        ["tmux", "display-message", "-p", "-t", target, "#{pane_dead}"],
        capture_output=True,
        text=True,
        check=False,
    )
    if res.returncode != 0:
        return False
    return res.stdout.strip() != "1"


def read_status(agent: str) -> dict | None:
    with closing(sqlite3.connect(DB_PATH)) as conn:
        row = conn.execute(
            "SELECT percent, activity, topic FROM agent_status WHERE agent=?",
            (agent,),
        ).fetchone()
    if not row:
        return None
    return {"percent": int(row[0]), "activity": row[1] or "", "topic": row[2] or ""}


def post_status(agent: str, percent: int, activity: str, topic: str) -> None:
    # heartbeat=True so room.py refreshes updated_at (sidebar freshness
    # tier) without advancing last_active_at — that field only moves on
    # real LLM/operator posts, which is what the chat-pane awake/asleep
    # badge keys off of.
    payload = json.dumps(
        {
            "agent": agent,
            "percent": percent,
            "activity": activity,
            "topic": topic,
            "heartbeat": True,
        }
    ).encode("utf-8")
    req = urllib.request.Request(
        ROOM_URL,
        data=payload,
        method="POST",
        headers={
            "Content-Type": "application/json",
            "Cookie": f"collab_room_token={urllib.parse.quote(ROOM_TOKEN)}",
        },
    )
    with urllib.request.urlopen(req, timeout=2.0) as resp:
        resp.read()


def tick(panes: list[tuple[str, str]]) -> None:
    for agent, target in panes:
        existing = read_status(agent)
        if pane_alive(target):
            activity = existing["activity"] if existing else "idle"
        else:
            activity = "(pane dead)"
        percent = existing["percent"] if existing else 0
        topic = existing["topic"] if existing else ""
        try:
            post_status(agent, percent, activity, topic)
        except Exception as exc:  # noqa: BLE001
            print(f"heartbeat post failed for {agent}: {exc}", file=sys.stderr, flush=True)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Heartbeat daemon for ACS Collab Workroom sidebar.")
    parser.add_argument(
        "--workroom",
        default=os.environ.get("ACS_WORKROOM", ""),
        help="Workroom name; agent IDs and tmux targets get prefixed. "
             "Default: $ACS_WORKROOM (empty = legacy unprefixed mode).",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)
    workroom = args.workroom.strip() or None
    panes = build_panes(workroom)
    print(
        f"acs-collab heartbeat daemon: workroom={workroom or '<legacy>'} "
        f"interval={INTERVAL_S}s url={ROOM_URL} panes={[a for a, _ in panes]}",
        file=sys.stderr,
        flush=True,
    )
    while True:
        try:
            tick(panes)
        except Exception as exc:  # noqa: BLE001
            print(f"heartbeat tick error: {exc}", file=sys.stderr, flush=True)
        time.sleep(INTERVAL_S)


if __name__ == "__main__":
    main()
