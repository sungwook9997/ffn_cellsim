"""Create the recommended ACS Workroom tmux session layout.

The layout separates always-idle chat panes from long-running work panes:

  claude-chat   receives PI workroom wrappers
  codex-chat    receives PI workroom wrappers
  claude-work   optional long-running Claude implementation work
  codex-work    optional long-running Codex review/implementation work
  win-ssh       optional SSH pane to the Windows A5000 workstation
  heartbeat     optional sidebar heartbeat daemon (refreshes agent_status
                rows for the four LLM panes every 5 s)

Existing sessions are never killed or renamed.
"""

from __future__ import annotations

import argparse
import pathlib
import shlex
import shutil
import subprocess
from typing import Iterable


DEFAULT_REPO = pathlib.Path(__file__).resolve().parents[2]

def _executable(name: str, fallback: str | None = None) -> str:
    path = shutil.which(name)
    if path is not None:
        return shlex.quote(path)
    if fallback is not None and pathlib.Path(fallback).exists():
        return shlex.quote(fallback)
    return shlex.quote(name)


CLAUDE_CMD = f"{_executable('claude')} --dangerously-skip-permissions"
CODEX_CMD = (
    f"{_executable('codex', '/Applications/Codex.app/Contents/Resources/codex')} "
    "--dangerously-bypass-approvals-and-sandbox"
)
WIN_SSH_CMD = "ssh win"
HEARTBEAT_CMD = "python3 tools/collab_mcp/heartbeat_daemon.py"
INTERACTIVE_SESSIONS = {
    "claude-chat",
    "codex-chat",
    "claude-work",
    "codex-work",
    "win-ssh",
}


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


def ensure_session(name: str, cwd: pathlib.Path, command: str | None = None) -> bool:
    if has_session(name):
        return False
    args = ["new-session", "-d", "-s", name, "-c", str(cwd)]
    if command is not None:
        args.append(command)
    tmux(*args)
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
    tmux("send-keys", "-t", target, "C-c")
    tmux("send-keys", "-t", target, line, "Enter")


def start_cli_if_shell(target: str, command: str) -> bool:
    current = pane_command(target)
    if current not in {"zsh", "bash", "fish", "sh"}:
        return False
    send_line(target, command)
    return True


def setup(
    cwd: pathlib.Path,
    start_chat: bool,
    start_work: bool,
    start_win_ssh: bool,
    start_heartbeat: bool,
) -> list[str]:
    actions: list[str] = []
    start_commands: dict[str, tuple[str, str]] = {}
    if start_chat:
        start_commands.update({
            "claude-chat": (CLAUDE_CMD, "started Claude in claude-chat"),
            "codex-chat": (CODEX_CMD, "started Codex in codex-chat"),
        })
    if start_work:
        start_commands.update({
            "claude-work": (CLAUDE_CMD, "started Claude in claude-work"),
            "codex-work": (CODEX_CMD, "started Codex in codex-work"),
        })
    if start_win_ssh:
        start_commands["win-ssh"] = (WIN_SSH_CMD, "started Windows SSH in win-ssh")
    if start_heartbeat:
        start_commands["heartbeat"] = (HEARTBEAT_CMD, "started heartbeat daemon")

    created_sessions: set[str] = set()
    sessions = [
        "claude-chat", "codex-chat", "claude-work", "codex-work",
        "win-ssh", "heartbeat",
    ]
    for session in sessions:
        command_and_action = start_commands.get(session)
        run_direct = command_and_action is not None and session not in INTERACTIVE_SESSIONS
        command = command_and_action[0] if run_direct else None
        created = ensure_session(session, cwd, command)
        actions.append(f"{'created' if created else 'exists'} {session}")
        if created:
            created_sessions.add(session)
            if run_direct and command_and_action is not None:
                actions.append(command_and_action[1])

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

    if start_win_ssh:
        if start_cli_if_shell("win-ssh:0.0", WIN_SSH_CMD):
            actions.append("started Windows SSH in win-ssh")

    if start_heartbeat:
        if "heartbeat" not in created_sessions and start_cli_if_shell("heartbeat:0.0", HEARTBEAT_CMD):
            actions.append("started heartbeat daemon")

    return actions


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Set up ACS Workroom tmux sessions.")
    parser.add_argument("--cwd", default=str(DEFAULT_REPO))
    parser.add_argument("--start-chat", action="store_true", help="Launch Claude/Codex in chat sessions.")
    parser.add_argument("--start-work", action="store_true", help="Launch Claude/Codex in work sessions.")
    parser.add_argument("--start-win-ssh", action="store_true", help="Launch ssh win in the Windows SSH session.")
    parser.add_argument("--start-heartbeat", action="store_true", help="Launch the sidebar heartbeat daemon.")
    return parser


def main(argv: Iterable[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)
    cwd = pathlib.Path(args.cwd).expanduser().resolve()
    actions = setup(
        cwd,
        start_chat=args.start_chat,
        start_work=args.start_work,
        start_win_ssh=args.start_win_ssh,
        start_heartbeat=args.start_heartbeat,
    )
    for action in actions:
        print(action)
    print("\nRelay targets:")
    print("  COLLAB_TMUX_CLAUDE_TARGET=claude-chat:0.0")
    print("  COLLAB_TMUX_CODEX_TARGET=codex-chat:0.0")
    print("  Windows SSH target: win-ssh:0.0")
    print("  Heartbeat target:   heartbeat:0.0")


if __name__ == "__main__":
    main()
