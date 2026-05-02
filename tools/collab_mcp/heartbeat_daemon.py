"""Heartbeat daemon for the ACS Collab Workroom sidebar.

Every 30 s, refresh the `agent_status` row for each of the four LLM panes
(`claude-chat`, `claude-work`, `codex-chat`, `codex-work`) so the sidebar
"Agents · 4 panes" panel stays visibly alive even when the pane is idle
between PI prompts. The daemon never invents progress numbers — it
re-POSTs the existing row's percent/activity unchanged, which only bumps
``updated_at`` and resets the freshness tick.

If the underlying tmux pane is dead, the daemon overwrites the row's
activity with ``"(pane dead)"`` and leaves the percent intact so PI can
spot a crashed CLI from the sidebar without inventing fake progress.
"""

from __future__ import annotations

import json
import os
import pathlib
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
INTERVAL_S = float(os.environ.get("COLLAB_HEARTBEAT_INTERVAL", "30"))

PANES = [
    ("claude-chat", "claude-chat:0.0"),
    ("claude-work", "claude-work:0.0"),
    ("codex-chat", "codex-chat:0.0"),
    ("codex-work", "codex-work:0.0"),
]


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


def tick() -> None:
    for agent, target in PANES:
        existing = read_status(agent)
        if not existing:
            continue
        if pane_alive(target):
            activity = existing["activity"]
        else:
            activity = "(pane dead)"
        try:
            post_status(agent, existing["percent"], activity, existing["topic"])
        except Exception as exc:  # noqa: BLE001
            print(f"heartbeat post failed for {agent}: {exc}", file=sys.stderr, flush=True)


def main() -> None:
    print(
        f"acs-collab heartbeat daemon: interval={INTERVAL_S}s url={ROOM_URL}",
        file=sys.stderr,
        flush=True,
    )
    while True:
        try:
            tick()
        except Exception as exc:  # noqa: BLE001
            print(f"heartbeat tick error: {exc}", file=sys.stderr, flush=True)
        time.sleep(INTERVAL_S)


if __name__ == "__main__":
    main()
