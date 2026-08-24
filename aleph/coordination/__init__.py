r"""Coordination primitives for several sessions sharing one repository and one GPU.

These are not physics and they own no simulation state.  They exist because the project moved from a
single Lead session to several concurrent ones, and the two things that do not survive that move are
(a) exclusive use of the single A5000 and (b) shared status files that every session wants to write.

The design principle throughout is **isolation by declared ownership, not by copying the tree**.  The
previous answer was a git worktree per parallel worker, which produced sibling checkouts (`..._b`) and
moved the merge problem rather than solving it.  Ownership declared in one file, enforced at commit
time, gives the same isolation with one working tree.

Modules:
    * :mod:`~aleph.coordination.gpu_lease` — a single-holder, self-expiring lease on the
      device, plus the append-only run index that answers "has this exact (build, config) already
      been measured?" before a session spends an hour re-measuring it.
    * :mod:`~aleph.coordination.ownership` — which session may commit which paths, and
      the check that enforces it.
"""

from __future__ import annotations

from aleph.coordination.gpu_lease import (
    DEFAULT_LEASE_PATH,
    DEFAULT_RUN_INDEX_PATH,
    Lease,
    LeaseHeld,
    RunIndex,
    RunRecord,
    acquire,
    read_lease,
    release,
)
from aleph.coordination.ownership import (
    OwnershipMap,
    OwnershipViolation,
    check_paths,
    load_ownership,
)

__all__ = [
    "DEFAULT_LEASE_PATH",
    "DEFAULT_RUN_INDEX_PATH",
    "Lease",
    "LeaseHeld",
    "OwnershipMap",
    "OwnershipViolation",
    "RunIndex",
    "RunRecord",
    "acquire",
    "check_paths",
    "load_ownership",
    "read_lease",
    "release",
]
