r"""Which session may commit which paths — isolation without a second working tree.

WHY.  Parallel work used to mean a git worktree per worker, which produced sibling checkouts beside
the repository and moved the integration problem instead of removing it.  The alternative is to keep
ONE tree and make the isolation a declaration: each concurrent session claims a set of path globs,
and a commit touching anything outside its claim is refused.

WHAT THE ENFORCEMENT IS FOR, precisely.  Not merge conflicts — git already reports those.  It is for
the silent case, which is the one that actually cost this project time: two sessions editing the same
shared status file minutes apart, where the loser's edit is not conflicted but simply absent, or is
swept into the winner's commit and attributed to it.  That happened on ``STATE.md`` twice in one hour
on 2026-07-28, and once a whole commit's staged content was carried away by a neighbouring session.

THE SHARED FILES ARE THE POINT, so they get their own category.  ``results_manifest.yaml`` and
``coverage_baseline.yaml`` are wanted by every session by construction: every run produces a claim, and
every claim belongs in both.  Declaring one owner for them would serialise the project on that session.
They are therefore marked ``shared`` — writable by nobody directly, and reached only through an
append-only source that renders into them.  A session that edits a shared file by hand is told which
append-only file to edit instead.

``STATE.md`` is deliberately NOT among them; it is ``free``.  Its tier-(a) block is GENERATED, and
``render_state_rows.py --check`` in the pre-commit hook already refuses a commit whose block drifted from
``aleph/docs/state_rows.yaml`` — so the part concurrent writers actually destroyed is protected
mechanically, while the prose sections are protected by discipline (read-then-write in one go).  A rule
refusing every ``STATE.md`` edit would have been refused-by-default rather than correct.  ``ownership.yaml``
carries the same reasoning at the ``free:`` key; if the two ever disagree, the YAML is the config that runs.

DELIBERATELY NOT DONE.  No locking, no ordering, no automatic merge.  A session that needs a path it
does not own asks for it; the failure mode of automatic transfer is that ownership means nothing.

Sanity Gate:
    * boundary: an empty claim list means the session owns NOTHING and every path is refused — a
      typo in a session name fails loudly rather than granting silent write access to the whole tree.
    * conservation/invariant: a path matched by two different sessions' globs is a configuration
      ERROR and is reported at load time, not at commit time, because overlapping ownership is
      exactly the state the file exists to prevent.
    * measurement-protocol: matching is on repo-relative POSIX paths, so the check gives the same
      answer regardless of where it is invoked from.
    * dimensional / sign-sense / numerical: not applicable — no quantities here.
"""

from __future__ import annotations

import fnmatch
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

__all__ = [
    "DEFAULT_OWNERSHIP_PATH",
    "OwnershipMap",
    "OwnershipViolation",
    "check_paths",
    "load_ownership",
]

#: Checked into the repo so the claim set is reviewable in a diff like anything else.
DEFAULT_OWNERSHIP_PATH = Path("ownership.yaml")


class OwnershipViolation(RuntimeError):
    """Raised when a session touches paths it has not claimed."""


@dataclass(frozen=True, slots=True)
class OwnershipMap:
    """The parsed claim set.

    Attributes:
        sessions: Session name -> list of repo-relative glob patterns it owns.
        shared: Globs no session may edit directly, mapped to the append-only file to edit instead.
        free: Globs any session may edit — scratch space and per-run output directories, which are
            unique per run by construction and so cannot collide.
    """

    sessions: dict[str, list[str]] = field(default_factory=dict)
    shared: dict[str, str] = field(default_factory=dict)
    free: list[str] = field(default_factory=list)

    def owner_of(self, path: str) -> str | None:
        """Return the session owning ``path``, or ``None`` if unclaimed.

        Args:
            path: Repo-relative POSIX path.

        Returns:
            The owning session name, or ``None``.
        """
        for session, globs in self.sessions.items():
            if any(fnmatch.fnmatch(path, g) for g in globs):
                return session
        return None

    def shared_redirect(self, path: str) -> str | None:
        """Return the append-only file to edit instead, if ``path`` is shared."""
        for pattern, redirect in self.shared.items():
            if fnmatch.fnmatch(path, pattern):
                return redirect
        return None

    def is_free(self, path: str) -> bool:
        """Whether ``path`` is in the unclaimed-by-design free space."""
        return any(fnmatch.fnmatch(path, g) for g in self.free)


def load_ownership(
    path: Path | str | None = None, *, data: Mapping[str, Any] | None = None
) -> OwnershipMap:
    """Load and validate the ownership declaration.

    Args:
        path: YAML file; defaults to :data:`DEFAULT_OWNERSHIP_PATH`.
        data: Pre-parsed mapping, for callers that already have it (and for tests, which then need
            no file at all).

    Returns:
        The validated :class:`OwnershipMap`.

    Raises:
        OwnershipViolation: If two sessions claim an identical glob — overlapping ownership defeats
            the purpose, so it is rejected at load rather than discovered at commit.
        FileNotFoundError: If the file is missing and no ``data`` was supplied.
    """
    if data is None:
        import yaml

        target = Path(path) if path is not None else DEFAULT_OWNERSHIP_PATH
        data = yaml.safe_load(target.read_text()) or {}

    sessions = {str(k): [str(g) for g in (v or [])] for k, v in (data.get("sessions") or {}).items()}
    seen: dict[str, str] = {}
    for session, globs in sessions.items():
        for g in globs:
            if g in seen:
                raise OwnershipViolation(
                    f"glob {g!r} is claimed by both {seen[g]!r} and {session!r}; overlapping "
                    "ownership is the state this file exists to prevent"
                )
            seen[g] = session
    return OwnershipMap(
        sessions=sessions,
        shared={str(k): str(v) for k, v in (data.get("shared") or {}).items()},
        free=[str(g) for g in (data.get("free") or [])],
    )


def check_paths(session: str, paths: list[str], ownership: OwnershipMap) -> list[str]:
    """Return one human-readable complaint per path ``session`` may not write.

    Args:
        session: The committing session's declared identity.
        paths: Repo-relative POSIX paths being committed.
        ownership: The loaded claim set.

    Returns:
        Complaints, empty when the commit is allowed. Returning rather than raising lets a caller
        report EVERY offending path at once — a hook that fails on the first one makes the author
        rediscover the rule file by file.
    """
    if session not in ownership.sessions:
        return [
            f"session {session!r} is not declared in the ownership file. Declare it (or unset "
            f"FFN_SESSION to opt out of parallel-session enforcement entirely); a session with no "
            f"claim owns nothing, which is the safe reading of a typo."
        ]
    complaints: list[str] = []
    for path in paths:
        if ownership.is_free(path):
            continue
        redirect = ownership.shared_redirect(path)
        if redirect is not None:
            complaints.append(
                f"{path} is SHARED — every session wants it, so no session edits it directly. "
                f"Edit {redirect} instead and re-render."
            )
            continue
        owner = ownership.owner_of(path)
        if owner is None:
            complaints.append(
                f"{path} is unclaimed. Add it to a session's claim in the ownership file, or to "
                f"`free:` if it is genuinely per-run output that cannot collide."
            )
        elif owner != session:
            complaints.append(f"{path} belongs to {owner!r}, not {session!r}.")
    return complaints
