"""Create the recommended ACS Workroom tmux session layout.

The active layout is one Claude + one Codex per workroom:

  <wr>-claude      receives workroom wrappers and may do the work directly
  <wr>-codex       receives workroom wrappers and may do the work directly
  <wr>-heartbeat   optional sidebar heartbeat daemon

`<wr>` is the workroom prefix (default `design-discussion`). The shared MCP
server (`mcp`) and room.py UI (`room`) are room-agnostic singletons and are
NOT created here — the parent launcher handles them.

Per-workroom MCP isolation: every Claude/Codex CLI is launched with the
workroom name baked into its MCP client config so the server-side
``X-Collab-Room`` header reaches the MCP server on every call. Without
this, the server-side ``room`` filter has nothing to filter on and both
workrooms write to ``design-discussion`` by default.

Existing sessions are never killed or renamed.
"""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import re
import shlex
import shutil
import subprocess
from typing import Iterable


DEFAULT_REPO = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_WORKROOM = os.environ.get("ACS_WORKROOM", "design-discussion")
WORKROOM_RE = re.compile(r"^[a-z0-9][a-z0-9-]{0,31}$")

def _executable(name: str, fallback: str | None = None) -> str:
    path = shutil.which(name)
    if path is not None:
        return shlex.quote(path)
    if fallback is not None and pathlib.Path(fallback).exists():
        return shlex.quote(fallback)
    return shlex.quote(name)


CLAUDE_BIN = _executable('claude', str(pathlib.Path.home() / '.local/bin/claude'))
CODEX_BIN = _executable('codex', '/Applications/Codex.app/Contents/Resources/codex')
WIN_SSH_CMD = "ssh win"
HEARTBEAT_CMD_TEMPLATE = "python3 tools/collab_mcp/heartbeat_daemon.py --workroom {workroom}"
DEFAULT_MCP_URL = os.environ.get("COLLAB_MCP_URL", "http://127.0.0.1:7878/mcp/")
PER_ROOM_LOG_DIR = pathlib.Path(os.environ.get("COLLAB_LOG_DIR_BASE", "/tmp/acs-collab"))
INTERACTIVE_LOGICAL = {"claude", "codex", "win-ssh"}


def write_claude_mcp_config(workroom: str, log_dir: pathlib.Path, token: str) -> pathlib.Path:
    """Write a per-workroom Claude MCP config and return its path.

    The launcher invokes Claude CLI with ``--mcp-config <path>
    --strict-mcp-config`` so the per-pane file is the only MCP source —
    no merge with the global ``~/.claude.json`` (which would otherwise
    overwrite the room header). The token is read from the launcher env
    at write time (``COLLAB_MCP_TOKEN``); the file is created mode 0600
    so a multi-user box does not leak the bearer.
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    cfg = {
        "mcpServers": {
            "acs-collab": {
                "type": "http",
                "url": DEFAULT_MCP_URL,
                "headers": {
                    "Authorization": f"Bearer {token}",
                    "X-Collab-Author": "claude",
                    "X-Collab-Room": workroom,
                },
            }
        }
    }
    path = log_dir / "claude_mcp.json"
    path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    try:
        path.chmod(0o600)
    except OSError:
        pass
    return path


def claude_command(workroom: str, mcp_config_path: pathlib.Path) -> str:
    """Return the shell line that launches Claude CLI in a workroom pane.

    Includes ``--model 'claude-opus-4-7[1m]' --effort xhigh`` to match
    desktop-app reasoning depth — the CLI default falls to a smaller
    effort tier and PI flagged the resulting "shallow" responses on
    2026-05-03 (memory: workroom_claude_flags).
    """
    return (
        f"ACS_WORKROOM={shlex.quote(workroom)} {CLAUDE_BIN} "
        f"--dangerously-skip-permissions "
        f"--model 'claude-opus-4-7[1m]' --effort xhigh "
        f"--mcp-config {shlex.quote(str(mcp_config_path))} --strict-mcp-config"
    )


def codex_command(workroom: str) -> str:
    """Return the shell line that launches Codex CLI in a workroom pane.

    Codex's global config supplies the URL and ``bearer_token_env_var``;
    we override only ``http_headers`` via ``-c`` so the X-Collab-Room
    header reaches the MCP server. The override replaces the entire
    ``http_headers`` table, so ``X-Collab-Author=codex`` is restated.
    """
    headers_toml = (
        f'mcp_servers.acs-collab.http_headers='
        f'{{ "X-Collab-Author" = "codex", "X-Collab-Room" = "{workroom}" }}'
    )
    return (
        f"ACS_WORKROOM={shlex.quote(workroom)} {CODEX_BIN} "
        f"--dangerously-bypass-approvals-and-sandbox "
        f"-c {shlex.quote(headers_toml)}"
    )


def validate_workroom(name: str) -> str:
    """Reject workroom names that would break tmux session naming or paths."""
    if not WORKROOM_RE.match(name):
        raise ValueError(
            f"invalid workroom name {name!r}: must match {WORKROOM_RE.pattern} "
            "(lowercase alphanum + dash, 1-32 chars, must start with alphanum)"
        )
    return name


def session_name(workroom: str, logical: str) -> str:
    """Map a logical pane (claude/codex/...) to its workroom-prefixed session."""
    return f"{workroom}-{logical}"


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
    result = tmux(*args, check=False)
    if result.returncode == 0:
        return True
    # Same-room launchers can race: both see no session, one creates it,
    # the other gets tmux's duplicate-session failure. Treat the loser as
    # idempotent success if the target now exists.
    if has_session(name):
        return False
    raise subprocess.CalledProcessError(
        result.returncode,
        ["tmux", *args],
        output=result.stdout,
        stderr=result.stderr,
    )


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
    workroom: str,
    start_chat: bool,
    start_work: bool,
    start_win_ssh: bool,
    start_heartbeat: bool,
) -> list[str]:
    workroom = validate_workroom(workroom)
    actions: list[str] = []
    heartbeat_cmd = HEARTBEAT_CMD_TEMPLATE.format(workroom=shlex.quote(workroom))
    log_dir = PER_ROOM_LOG_DIR / workroom
    token = os.environ.get("COLLAB_MCP_TOKEN", "")
    start_agents = start_chat or start_work
    if start_agents and not token:
        # Without a token we cannot generate a usable per-workroom Claude
        # MCP config. Fail loud rather than write an unauthenticated
        # config that would silently break the room header propagation.
        raise RuntimeError(
            "COLLAB_MCP_TOKEN is unset; cannot write per-workroom Claude MCP "
            "config. The launcher should export it before invoking "
            "setup_tmux_workroom."
        )
    claude_cfg = (
        write_claude_mcp_config(workroom, log_dir, token)
        if start_agents and token
        else None
    )
    claude_cmd = claude_command(workroom, claude_cfg) if claude_cfg else CLAUDE_BIN
    codex_cmd = codex_command(workroom)
    start_commands: dict[str, tuple[str, str]] = {}
    if start_agents:
        start_commands.update({
            "claude": (claude_cmd, f"started Claude in {workroom}-claude"),
            "codex": (codex_cmd, f"started Codex in {workroom}-codex"),
        })
    if start_win_ssh:
        start_commands["win-ssh"] = (WIN_SSH_CMD, f"started Windows SSH in {workroom}-win-ssh")
    if start_heartbeat:
        start_commands["heartbeat"] = (heartbeat_cmd, f"started heartbeat daemon ({workroom})")

    created_sessions: set[str] = set()
    logical_sessions = ["claude", "codex"]
    if start_win_ssh:
        logical_sessions.append("win-ssh")
    if start_heartbeat:
        logical_sessions.append("heartbeat")
    for logical in logical_sessions:
        physical = session_name(workroom, logical)
        command_and_action = start_commands.get(logical)
        run_direct = (
            command_and_action is not None
            and logical not in INTERACTIVE_LOGICAL
        )
        command = command_and_action[0] if run_direct else None
        created = ensure_session(physical, cwd, command)
        actions.append(f"{'created' if created else 'exists'} {physical}")
        if created:
            created_sessions.add(logical)
            if run_direct and command_and_action is not None:
                actions.append(command_and_action[1])

    if start_agents:
        if start_cli_if_shell(f"{session_name(workroom, 'claude')}:0.0", claude_cmd):
            actions.append(f"started Claude in {workroom}-claude")
        if start_cli_if_shell(f"{session_name(workroom, 'codex')}:0.0", codex_cmd):
            actions.append(f"started Codex in {workroom}-codex")

    if start_win_ssh:
        if start_cli_if_shell(f"{session_name(workroom, 'win-ssh')}:0.0", WIN_SSH_CMD):
            actions.append(f"started Windows SSH in {workroom}-win-ssh")

    if start_heartbeat:
        if (
            "heartbeat" not in created_sessions
            and start_cli_if_shell(f"{session_name(workroom, 'heartbeat')}:0.0", heartbeat_cmd)
        ):
            actions.append(f"started heartbeat daemon ({workroom})")

    return actions


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Set up ACS Workroom tmux sessions.")
    parser.add_argument("--cwd", default=str(DEFAULT_REPO))
    parser.add_argument(
        "--workroom",
        default=DEFAULT_WORKROOM,
        help="Workroom name; tmux sessions are prefixed with this. "
             "Default: $ACS_WORKROOM or 'design-discussion'.",
    )
    parser.add_argument("--start-chat", action="store_true", help="Launch Claude/Codex agent sessions.")
    parser.add_argument("--start-work", action="store_true", help="Deprecated alias for --start-chat.")
    parser.add_argument("--start-win-ssh", action="store_true", help="Launch ssh win in the Windows SSH session.")
    parser.add_argument("--start-heartbeat", action="store_true", help="Launch the sidebar heartbeat daemon.")
    return parser


def main(argv: Iterable[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)
    cwd = pathlib.Path(args.cwd).expanduser().resolve()
    workroom = validate_workroom(args.workroom)
    actions = setup(
        cwd,
        workroom=workroom,
        start_chat=args.start_chat,
        start_work=args.start_work,
        start_win_ssh=args.start_win_ssh,
        start_heartbeat=args.start_heartbeat,
    )
    for action in actions:
        print(action)
    print(f"\nWorkroom: {workroom}")
    print("Relay targets:")
    print(f"  COLLAB_TMUX_CLAUDE_TARGET={workroom}-claude:0.0")
    print(f"  COLLAB_TMUX_CODEX_TARGET={workroom}-codex:0.0")
    print(f"  Windows SSH target: {workroom}-win-ssh:0.0")
    print(f"  Heartbeat target:   {workroom}-heartbeat:0.0")


if __name__ == "__main__":
    main()
