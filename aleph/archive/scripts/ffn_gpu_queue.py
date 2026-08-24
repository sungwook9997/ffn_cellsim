#!/usr/bin/env python
r"""Keep the device fed: any lane appends a job, one drainer runs them back to back.

WHY THIS EXISTS.  The lease made the device safe to share and did nothing to keep it BUSY.  Measured on
2026-07-28: jobs on this device range from 7 s (the T10 exterior-medium gate) to hours (a native resting
solve), and the only way a lane could react to contention was to give up — two lane reports that day
carry the same sentence, "native run 0, the A5000 was leased all session", while the holder sat at 5% of
device memory.  `--wait` fixed the giving-up.  This fixes the other half: between one job finishing and a
human noticing, the device is idle, and at these durations that gap is routinely longer than the jobs.

    # a lane, at any time, from anywhere in the repo
    python aleph/scripts/ffn_gpu_queue.py add --lane sf-motor --minutes 20 \
        --reason "sf_motor native gate at the new step scale" \
        --population "904 nodes / 2,712 DOF — SLICE" \
        -- aleph/scripts/ac_gate_b_sf_motor_native.py --k-axial 1000

    python aleph/scripts/ffn_gpu_queue.py list
    python aleph/scripts/ffn_gpu_queue.py drain          # runs everything, then exits
    python aleph/scripts/ffn_gpu_queue.py drain --follow # ...and keeps waiting for new jobs

ONE DRAINER, MANY PRODUCERS.  Appends are atomic single-line `O_APPEND` writes, which the kernel does not
interleave below `PIPE_BUF`, so lanes never coordinate with each other to enqueue.  Draining is guarded by
its own lock file so two drainers cannot both claim the same job; the lock carries a pid and expires, for
the same reason the device lease does.

WHAT THIS DELIBERATELY DOES NOT DO.  It does not bypass the lease, the duplicate-run index, or the remote
digest check — every job is launched through `ffn_gpu.py run`, so a queued job is subject to exactly the
same refusals as a hand-typed one, including "this (build, config) already ran".  A queue that could skip
those would turn one wrong record into a batch of them.

It also does not reorder, prioritise, or pack jobs concurrently.  Ordering is arrival order because any
other policy is a scheduling decision nobody has made, and concurrency needs a per-run VRAM budget plus a
lease that can be shared — a coordination-contract change (`read_lease` fails OPEN, so changing the lease
schema while sessions are live makes every running session see FREE at once).

Sanity Gate:
    * boundary: an empty queue drains to a no-op; a job whose driver does not resolve is refused by
      `ffn_gpu.py run` before the lease is taken, so a typo cannot idle the device.
    * conservation/invariant: a job is moved to the done-log only after its runner exits, and the exit
      code is recorded with it, so a crashed drainer leaves the job pending rather than silently dropped.
    * measurement-protocol: each job records its own lane, reason and population, and those are passed
      through to the run record — a queued run is not anonymous.
    * numerical / dimensional / sign-sense: not applicable — this module computes no quantity.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

#: Queue of pending jobs, one JSON object per line, appended by any lane.
DEFAULT_QUEUE = Path.home() / ".ffn" / "queue.jsonl"
#: Append-only record of what ran, with exit codes.
DEFAULT_DONE = Path.home() / ".ffn" / "queue_done.jsonl"
#: Held while a drainer is running so two do not claim the same job.
DEFAULT_LOCK = Path.home() / ".ffn" / "queue.lock"

#: How long a drainer's lock stays valid without refresh. Long enough to cover the longest job seen
#: (a native resting solve), short enough that a dead drainer does not wedge the queue for a day.
LOCK_TTL_SECONDS = 6 * 3600.0
#: Poll interval when following an empty queue.
FOLLOW_POLL_SECONDS = 15.0


def _session() -> str:
    return os.environ.get("FFN_SESSION") or f"unnamed@pid{os.getpid()}"


def append_job(job: dict, path: Path = DEFAULT_QUEUE) -> None:
    """Append one job atomically.

    A single ``O_APPEND`` write below ``PIPE_BUF`` is not interleaved by the kernel, so concurrent
    lanes need no lock to enqueue — which is the point, since a lane that had to coordinate to add
    work would be back to coordinating for the device.

    Args:
        job: The job record; must be JSON-serialisable and fit on one line.
        path: Queue file.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(job, separators=(",", ":")) + "\n"
    with open(path, "a", encoding="utf-8") as handle:
        handle.write(line)


def read_jobs(path: Path = DEFAULT_QUEUE) -> list[dict]:
    """Return pending jobs in arrival order, skipping any line that will not parse.

    A malformed line is skipped rather than fatal: one bad append must not block every other lane's
    work, and the line stays in the file for a human to look at.
    """
    if not path.exists():
        return []
    jobs = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not raw.strip():
            continue
        try:
            jobs.append(json.loads(raw))
        except ValueError:
            continue
    return jobs


def _drainer_alive(pid: int) -> bool:
    """Whether the recorded drainer process still exists on THIS machine."""
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)          # signal 0 tests existence without touching the process
    except ProcessLookupError:
        return False
    except PermissionError:
        return True              # exists, owned by someone else
    return True


def _take_lock(path: Path = DEFAULT_LOCK) -> bool:
    """Claim the drainer role; return False only if a drainer is REALLY still running.

    Liveness is decided by the recorded pid, with the TTL as a backstop — the same order `gpu_lease`
    uses, and for the reason it learned: a drainer killed with SIGKILL never runs its `finally`, and a
    TTL-only check then wedges the queue for the whole TTL. Measured on 2026-07-28: killing a drainer
    left three native runs queued and unstartable behind a six-hour lock, and the ETA reported for them
    was wrong because nothing had begun.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    now = time.time()
    try:
        payload = json.loads(path.read_text())
        pid = int(payload.get("pid", 0))
        fresh = now - float(payload.get("at", 0.0)) < LOCK_TTL_SECONDS
        if fresh and _drainer_alive(pid):
            return False
        if fresh:
            print(f"[queue] taking over from drainer pid {pid}, which is gone", file=sys.stderr)
    except (OSError, ValueError, TypeError):
        pass
    path.write_text(json.dumps({"holder": _session(), "pid": os.getpid(), "at": now}) + "\n")
    return True


def _release_lock(path: Path = DEFAULT_LOCK) -> None:
    try:
        payload = json.loads(path.read_text())
        if payload.get("pid") == os.getpid():
            path.unlink()
    except (OSError, ValueError, TypeError):
        pass


#: Exit codes from `ffn_gpu.py run` that mean the job NEVER STARTED and the cause is repairable:
#: 3 = this (build, config) already ran, 4 = the remote tree does not hold the code the run would stamp.
#: A job refused for either reason must go back on the queue, not into the done log — it has no result to
#: record, and the drainer is usually unattended, so discarding it loses the work silently.
RETRYABLE_REFUSALS: frozenset[int] = frozenset({3, 4})

#: Bounded so a permanently-broken job cannot spin the queue forever; after this it is recorded as done
#: with its refusal code, which is a visible failure rather than an invisible one.
MAX_REQUEUES: int = 3

#: Enough for a human to rsync between retries, short enough not to idle the device on a transient.
REQUEUE_BACKOFF_SECONDS: float = 120.0


def _pop_first(path: Path = DEFAULT_QUEUE) -> dict | None:
    """Remove and return the first pending job, or None."""
    jobs = read_jobs(path)
    if not jobs:
        return None
    head, rest = jobs[0], jobs[1:]
    tmp = path.with_suffix(f".{os.getpid()}.tmp")
    tmp.write_text("".join(json.dumps(j, separators=(",", ":")) + "\n" for j in rest))
    os.replace(tmp, path)
    return head


def run_job(job: dict, repo: Path, *, wait_minutes: float) -> int:
    """Launch one job through `ffn_gpu.py run`, inheriting every guard that path applies."""
    argv = [sys.executable, str(repo / "aleph" / "scripts" / "ffn_gpu.py"), "run",
            "--minutes", str(job.get("minutes", 30)),
            "--wait", str(wait_minutes),
            "--reason", job.get("reason", "") or f"queued by {job.get('lane', '?')}"]
    if job.get("population"):
        argv += ["--population", job["population"]]
    if job.get("artifact"):
        argv += ["--artifact", job["artifact"]]
    if job.get("force"):
        argv += ["--force"]
    argv += ["--"] + list(job.get("command", []))

    env = dict(os.environ)
    if job.get("lane"):
        env["FFN_SESSION"] = job["lane"]
    print(f"[queue] -> {job.get('lane', '?')}: {' '.join(job.get('command', [])[:2])}", flush=True)
    return subprocess.run(argv, cwd=repo, env=env, check=False).returncode


def reject_malformed_argv(command: list[str], repo: Path) -> str:
    """Return why this argv is structurally wrong, or ``""`` if it looks like a real argv.

    WHY IT IS WORTH CHECKING AT ENQUEUE.  A malformed entry costs a GPU SLOT to discover: the drainer takes
    the lease, ssh's to the run host, the driver's argparse exits 2, and the job is recorded as done.  That
    happened on 2026-07-29 because **zsh does not word-split an unquoted variable** — a whole flag string
    arrived as ONE argv element, ~200 characters long, and the driver saw a single unrecognised argument.

    WHAT IT CANNOT CATCH, said plainly so nobody trusts it further than it goes: a command that parses but is
    WRONG.  The same night, a re-queue typed from memory silently dropped ``--out``, ``--membrane-subdiv``
    and ``--k-xb``; every remaining token was a valid flag, the run executed, and it wrote nothing.  Only the
    post-run artifact probe in `ffn_gpu.py` catches that one.  These two guards cover different halves.
    """
    if not (repo / command[0]).is_file():
        return f"{command[0]} does not exist in this tree"
    for i, token in enumerate(command):
        if " " in token and token.lstrip().startswith("-"):
            return (f"argv[{i}] is a single {len(token)}-character token containing spaces and starting with "
                    f"'-': {token[:60]!r}...  That is one flag string that was never word-split — in zsh an "
                    f"unquoted $VAR does NOT split. Pass the flags literally, or use ${{=VAR}}")
    return ""


def cmd_add(args: argparse.Namespace) -> int:
    command = args.command[1:] if args.command and args.command[0] == "--" else args.command
    if not command:
        print("nothing to queue", file=sys.stderr)
        return 2
    malformed = reject_malformed_argv(command, Path(__file__).resolve().parents[2])
    if malformed:
        print(f"REFUSED — {malformed}", file=sys.stderr)
        return 2
    append_job({"lane": args.lane or _session(), "minutes": args.minutes, "reason": args.reason,
                "population": args.population, "artifact": args.artifact, "force": args.force,
                "command": command, "queued_at": time.time()})
    print(f"queued for {args.lane or _session()!r}: {' '.join(command[:3])}")
    return 0


def cmd_list(args: argparse.Namespace) -> int:
    jobs = read_jobs()
    if not jobs:
        print("queue is empty")
        return 0
    print(f"{len(jobs)} pending:")
    for i, job in enumerate(jobs, 1):
        print(f"  {i}. [{job.get('lane', '?')}] {' '.join(job.get('command', [])[:2])}"
              f"   ({job.get('minutes', '?')} min)  {job.get('reason', '')[:60]}")
    return 0


def cmd_drain(args: argparse.Namespace) -> int:
    repo = Path(__file__).resolve().parents[2]
    if not _take_lock():
        print("a drainer already holds the queue (lock is live) — not starting a second one",
              file=sys.stderr)
        return 2
    ran = 0
    try:
        while True:
            job = _pop_first()
            if job is None:
                if not args.follow:
                    break
                time.sleep(FOLLOW_POLL_SECONDS)
                continue
            code = run_job(job, repo, wait_minutes=args.wait)
            if code in RETRYABLE_REFUSALS and int(job.get("requeues", 0)) < MAX_REQUEUES:
                # The job never RAN — `ffn_gpu.py run` refused it before taking the device, and the
                # reason is repairable (a stale remote tree, or a duplicate (build, config)). Discarding
                # it here is silent data loss: on 2026-07-29 five sweep points were consumed this way
                # because source was committed between two drains, and the queue reported "drained"
                # while four of the six points a sweep needed no longer existed anywhere.
                job = {**job, "requeues": int(job.get("requeues", 0)) + 1,
                       "last_refusal_exit": code, "last_refusal_at": time.time()}
                append_job(job)
                print(f"[queue] <- exit {code} REFUSED (never ran) — requeued at the tail "
                      f"({job['requeues']}/{MAX_REQUEUES}). Fix the cause (usually: rsync source to the "
                      f"run host) and it will be retried.", flush=True)
                time.sleep(REQUEUE_BACKOFF_SECONDS)
                continue
            DEFAULT_DONE.parent.mkdir(parents=True, exist_ok=True)
            with open(DEFAULT_DONE, "a", encoding="utf-8") as handle:
                handle.write(json.dumps({**job, "exit_code": code, "finished_at": time.time()},
                                        separators=(",", ":")) + "\n")
            ran += 1
            if code in RETRYABLE_REFUSALS:
                print(f"[queue] <- exit {code} REFUSED and RETRIED {MAX_REQUEUES}x — giving up and "
                      f"recording it. This job is now gone from the queue.", flush=True)
            else:
                print(f"[queue] <- exit {code}", flush=True)
    finally:
        _release_lock()
    print(f"[queue] drained {ran} job(s)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("add", help="append a job; any lane, any time")
    p.add_argument("--lane", default="")
    p.add_argument("--minutes", type=float, default=30.0, help="lease length this job needs")
    p.add_argument("--reason", default="")
    p.add_argument("--population", default="", help="census string — a queued run is not anonymous")
    p.add_argument("--artifact", default="")
    p.add_argument("--force", action="store_true", help="run even if this (build, config) ran before")
    p.add_argument("command", nargs=argparse.REMAINDER)
    p.set_defaults(func=cmd_add)

    p = sub.add_parser("list", help="what is pending")
    p.set_defaults(func=cmd_list)

    p = sub.add_parser("drain", help="run pending jobs back to back")
    p.add_argument("--wait", type=float, default=120.0,
                   help="minutes each job may wait for the device before giving up")
    p.add_argument("--follow", action="store_true",
                   help="keep running after the queue empties, picking up jobs as lanes add them")
    p.set_defaults(func=cmd_drain)

    args = ap.parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
