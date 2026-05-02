"""Create the recommended ACS Workroom tmux session layout.

The layout separates always-idle chat panes from long-running work panes:

  claude-chat   receives PI workroom wrappers
  codex-chat    receives PI workroom wrappers
  claude-work   optional long-running Claude implementation work
  codex-work    optional long-running Codex review/implementation work

Existing sessions are never killed or renamed.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
from typing import Iterable


DEFAULT_REPO = pathlib.Path(__file__).resolve().parents[2]

CLAUDE_CMD = "claude --dangerously-skip-permissions"
CODEX_CMD = "codex --dangerously-bypass-approvals-and-sandbox"


def tmux(*args: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["tmux", *args],
        check=check,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def has_session(name: str) -> bool:
    return tmux("has-session", "-t", name, check=False).returncode == 0


def ensure_session(name: str, cwd: pathlib.Path) -> bool:
    if has_session(name):
        return False
    tmux("new-session", "-d", "-s", name, "-c", str(cwd))
    return True


def pane_command(target: str) -> str:
    result = tmux(
        "list-panes",
        "-t",
        target,
        "-F",
        "#{pane_current_command}",
    )
    return result.stdout.strip().splitlines()[0]


def send_line(target: str, line: str) -> None:
    tmux("send-keys", "-t", target, "-l", line)
    tmux("send-keys", "-t", target, "C-m")


def start_cli_if_shell(target: str, command: str) -> bool:
    current = pane_command(target)
    if current not in {"zsh", "bash", "fish", "sh"}:
        return False
    send_line(target, command)
    return True


def setup(cwd: pathlib.Path, start_chat: bool, start_work: bool) -> list[str]:
    actions: list[str] = []
    sessions = ["claude-chat", "codex-chat", "claude-work", "codex-work"]
    for session in sessions:
        created = ensure_session(session, cwd)
        actions.append(f"{'created' if created else 'exists'} {session}")

    if start_chat:
        if start_cli_if_shell("claude-chat:0.0", CLAUDE_CMD):
            actions.append("started Claude in claude-chat")
        if start_cli_if_shell("codex-chat:0.0", CODEX_CMD):
            actions.append("started Codex in codex-chat")

    if start_work:
        if start_cli_if_shell("claude-work:0.0", CLAUDE_CMD):
            actions.append("started Claude in claude-work")
        if start_cli_if_shell("codex-work:0.0", CODEX_CMD):
            actions.append("started Codex in codex-work")

    return actions


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Set up ACS Workroom tmux sessions.")
    parser.add_argument("--cwd", default=str(DEFAULT_REPO))
    parser.add_argument("--start-chat", action="store_true", help="Launch Claude/Codex in chat sessions.")
    parser.add_argument("--start-work", action="store_true", help="Launch Claude/Codex in work sessions.")
    return parser


def main(argv: Iterable[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)
    cwd = pathlib.Path(args.cwd).expanduser().resolve()
    actions = setup(cwd, start_chat=args.start_chat, start_work=args.start_work)
    for action in actions:
        print(action)
    print("\nRelay targets:")
    print("  COLLAB_TMUX_CLAUDE_TARGET=claude-chat:0.0")
    print("  COLLAB_TMUX_CODEX_TARGET=codex-chat:0.0")


if __name__ == "__main__":
    main()
