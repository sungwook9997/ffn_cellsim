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
    mkdir -p "$TARGET_DIR" /tmp/acs-collab
    sed "s|__REPO_ROOT__|$REPO_ROOT|g" "$TEMPLATE" > "$TARGET"
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
