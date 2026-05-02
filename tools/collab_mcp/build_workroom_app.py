"""Build a macOS .app bundle that wraps the ACS Workroom browser surface.

Stdlib-only, no Xcode needed. The bundle's executable is a thin zsh shim that
delegates to `tools/collab_mcp/launch_workroom_mac.command`, which is the
single source of truth for full-stack boot (syncthing check, tmux 4 work
sessions, relay, room.py, Chrome). Re-run this script only when the repo
moves — the .command file itself can be edited freely without rebuilding.

Default output: ~/Applications/ACSWorkroom.app — user-writable, no sudo.
"""

from __future__ import annotations

import argparse
import os
import plistlib
import shutil
import stat
import sys
from pathlib import Path

DEFAULT_APP_NAME = "ACSWorkroom"
DEFAULT_BUNDLE_ID = "com.activecellsim.workroom"
DEFAULT_OUTPUT = Path.home() / "Applications"

LAUNCHER_TEMPLATE = """#!/bin/zsh
set -euo pipefail

REPO_ROOT={repo_root_quoted}
exec "$REPO_ROOT/tools/collab_mcp/launch_workroom_mac.command"
"""


def _shell_quote(path: str) -> str:
    return "'" + path.replace("'", "'\\''") + "'"


def _detect_repo_root(script_path: Path) -> Path:
    return script_path.resolve().parents[2]


def _write_executable(path: Path, body: str) -> None:
    path.write_text(body, encoding="utf-8")
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _write_plist(path: Path, app_name: str, bundle_id: str) -> None:
    plist = {
        "CFBundleName": app_name,
        "CFBundleDisplayName": app_name,
        "CFBundleExecutable": app_name,
        "CFBundleIdentifier": bundle_id,
        "CFBundleVersion": "1.0",
        "CFBundleShortVersionString": "1.0",
        "CFBundlePackageType": "APPL",
        "CFBundleInfoDictionaryVersion": "6.0",
        "LSMinimumSystemVersion": "11.0",
        "NSHighResolutionCapable": True,
        "LSUIElement": False,
    }
    with path.open("wb") as f:
        plistlib.dump(plist, f)


def build(output_dir: Path, app_name: str, bundle_id: str, repo_root: Path) -> Path:
    if not repo_root.is_dir():
        raise SystemExit(f"repo root not a directory: {repo_root}")
    room_py = repo_root / "tools" / "collab_mcp" / "room.py"
    if not room_py.is_file():
        raise SystemExit(f"room.py not found at {room_py}")

    output_dir.mkdir(parents=True, exist_ok=True)
    app_path = output_dir / f"{app_name}.app"
    contents = app_path / "Contents"
    macos = contents / "MacOS"
    resources = contents / "Resources"

    if app_path.exists():
        shutil.rmtree(app_path)
    macos.mkdir(parents=True)
    resources.mkdir(parents=True)

    _write_plist(contents / "Info.plist", app_name, bundle_id)
    launcher_body = LAUNCHER_TEMPLATE.format(repo_root_quoted=_shell_quote(str(repo_root)))
    _write_executable(macos / app_name, launcher_body)

    return app_path


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT,
                        help=f"Directory to place the .app (default: {DEFAULT_OUTPUT})")
    parser.add_argument("--name", default=DEFAULT_APP_NAME,
                        help=f"App display name (default: {DEFAULT_APP_NAME})")
    parser.add_argument("--bundle-id", default=DEFAULT_BUNDLE_ID,
                        help=f"CFBundleIdentifier (default: {DEFAULT_BUNDLE_ID})")
    parser.add_argument("--repo-root", type=Path, default=None,
                        help="Override repo root (default: auto-detect from this script's location)")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv if argv is not None else sys.argv[1:])
    repo_root = args.repo_root or _detect_repo_root(Path(__file__))
    app_path = build(args.output.expanduser(), args.name, args.bundle_id, repo_root.expanduser())
    print(f"built: {app_path}")
    print(f"repo_root baked in: {repo_root}")
    print(f"launch: open {app_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
