#!/bin/zsh
# Install / uninstall the ACS Collab Workroom user LaunchAgent.
#
# Usage:
#   tools/collab_mcp/install_launchd.sh install    # default
#   tools/collab_mcp/install_launchd.sh uninstall
#   tools/collab_mcp/install_launchd.sh status
#
# The installer renders the plist template (substituting __REPO_ROOT__),
# drops it into ~/Library/LaunchAgents/, then `launchctl bootstrap`s it
# under the current GUI user so the workroom comes up on every login.

set -euo pipefail

SCRIPT_DIR="${0:A:h}"
REPO_ROOT="${SCRIPT_DIR:h:h}"
LABEL="com.activecellsim.workroom"
TEMPLATE="$SCRIPT_DIR/$LABEL.plist"
TARGET_DIR="$HOME/Library/LaunchAgents"
TARGET="$TARGET_DIR/$LABEL.plist"

action="${1:-install}"

case "$action" in
  install)
    if [[ ! -f "$TEMPLATE" ]]; then
      echo "template not found: $TEMPLATE" >&2
      exit 1
    fi
    if [[ ! -x "$REPO_ROOT/tools/collab_mcp/launch_workroom_mac.command" ]]; then
      echo "launcher not executable: $REPO_ROOT/tools/collab_mcp/launch_workroom_mac.command" >&2
      exit 1
    fi
    if [[ -z "${COLLAB_MCP_TOKEN:-}" ]]; then
      echo "COLLAB_MCP_TOKEN unset in current shell — refuse to install a plist with an empty token." >&2
      echo "Source ~/.zshenv (or export COLLAB_MCP_TOKEN=...) and re-run this installer." >&2
      exit 1
    fi
    if [[ ! -x "$REPO_ROOT/.venv-collab/bin/python" ]]; then
      echo "venv missing: $REPO_ROOT/.venv-collab/bin/python — bootstrap the collab venv before installing." >&2
      exit 1
    fi
    if ! "$REPO_ROOT/.venv-collab/bin/python" -c 'import mcp' >/dev/null 2>&1; then
      echo "$REPO_ROOT/.venv-collab/bin/python cannot import 'mcp' — run 'pip install mcp[cli]>=1.2' inside the venv before installing." >&2
      exit 1
    fi
    mkdir -p "$TARGET_DIR" /tmp/acs-collab
    # Template substitution keeps machine-specific paths out of the repo.
    # Token is read from the current shell environment, never written
    # back into the in-repo template, so the secret only lives in the
    # rendered file under ~/Library/LaunchAgents/.
    umask 077
    sed -e "s|__REPO_ROOT__|$REPO_ROOT|g" \
        -e "s|__HOME__|$HOME|g" \
        -e "s|__COLLAB_MCP_TOKEN__|$COLLAB_MCP_TOKEN|g" \
        "$TEMPLATE" > "$TARGET"
    chmod 600 "$TARGET"
    plutil -lint "$TARGET" >/dev/null

    uid="$(id -u)"
    domain="gui/$uid"
    # Idempotent: bootout if already loaded, then bootstrap.
    launchctl bootout "$domain" "$TARGET" 2>/dev/null || true
    launchctl bootstrap "$domain" "$TARGET"
    launchctl enable "$domain/$LABEL"

    echo "installed: $TARGET"
    echo "label:     $LABEL"
    echo "domain:    $domain"
    echo "logs:      /tmp/acs-collab/launchd.{out,err}"
    echo
    echo "Run on demand:    launchctl kickstart -k $domain/$LABEL"
    echo "Disable on login: launchctl disable $domain/$LABEL"
    echo "Uninstall:        $0 uninstall"
    ;;

  uninstall)
    uid="$(id -u)"
    domain="gui/$uid"
    if [[ -f "$TARGET" ]]; then
      launchctl bootout "$domain" "$TARGET" 2>/dev/null || true
      rm -f "$TARGET"
      echo "removed: $TARGET"
    else
      echo "not installed (no $TARGET)"
    fi
    ;;

  status)
    uid="$(id -u)"
    domain="gui/$uid"
    if [[ -f "$TARGET" ]]; then
      echo "plist:  $TARGET"
    else
      echo "plist:  (not installed)"
    fi
    if launchctl print "$domain/$LABEL" >/dev/null 2>&1; then
      echo "state:  loaded"
      launchctl print "$domain/$LABEL" | grep -E '^\s*(state|last exit code|program arguments)' || true
    else
      echo "state:  not loaded"
    fi
    ;;

  *)
    echo "usage: $0 {install|uninstall|status}" >&2
    exit 2
    ;;
esac
