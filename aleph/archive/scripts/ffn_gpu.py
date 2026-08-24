#!/usr/bin/env python
r"""One CLI for the shared GPU: take the lease, refuse a duplicate run, record what executed.

There is one A5000 and several concurrent sessions. Before this, the protocol was "whoever launches
first wins" — two native runs would land on the same 16 GB device, and neither artifact would record
that it had shared it. A run that was throttled or swapped by a neighbour is not the run its config
claims to be, and nothing downstream could tell.

    # who has it?
    python aleph/scripts/ffn_gpu.py status

    # take it for the length of the job, saying what for
    python aleph/scripts/ffn_gpu.py acquire --minutes 45 --reason "sf_motor native gate"

    # run natively WITH the lease enforced, the duplicate check applied, and the run recorded
    python aleph/scripts/ffn_gpu.py run --minutes 45 --reason "sf_motor native gate" -- \
        aleph/scripts/ac_gate_b_sf_motor_native.py --k-axial 1000 ...

    python aleph/scripts/ffn_gpu.py release
    python aleph/scripts/ffn_gpu.py index --limit 20

WHY THE LEASE LIVES ON THIS MACHINE. Every session that drives the GPU runs here and reaches gbook
over ssh, so the contenders are local even though the device is not. A local lease needs no network
round-trip to be correct and cannot be defeated by ssh being down. It does NOT cover someone logged
into gbook running a job by hand — that case is a conversation, not a lock.

WHY `run` REFUSES A DUPLICATE. The engine's run records already stamp the build commit and a
canonical config hash, so "has this exact configuration already been measured on this exact source?"
has been mechanically answerable for some time and nothing asked it. The observed consequence was
sessions re-running work that existed and overwriting each other's output. `run` asks before
spending the device, and `--force` records the duplicate as deliberate rather than hiding it.

Sanity Gate:
    * boundary: a run whose command does not resolve is rejected before the lease is taken, so a typo
      does not idle the device for the lease duration.
    * conservation/invariant: the lease is released in a `finally`, so a failing run frees the device;
      the expiry is the backstop for the case where the process dies without unwinding.
    * measurement-protocol: the duplicate check is on (build commit, config hash) TOGETHER — the same
      config on a different build is a different run, and the index records the population string so a
      reader sees at a glance whether a result came from native or from a slice.
    * numerical / dimensional / sign-sense: not applicable — this module computes no quantity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

# This CLI is meant to be run directly (`python aleph/scripts/ffn_gpu.py status`), and outside
# pytest there is no path setup, so importing the engine would fail. Prepending the repo root is the
# minimum that makes it work from anywhere; note it imports `aleph.…`, NOT the bare `ac.…` identity
# that 22 other scripts still use and that double-loads modules outside pytest (STATE.md (a)).
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from aleph.coordination.gpu_lease import (  # noqa: E402
    LeaseHeld,
    RunIndex,
    RunRecord,
    acquire,
    read_lease,
    release,
)

#: How a session names itself. Falls back to something a human can still trace.
SESSION_ENV = "FFN_SESSION"

#: Poll interval while ``--wait`` blocks. Short enough that a 7-second job like the T10 gate is not
#: left queueing long behind its predecessor, long enough not to spin on the filesystem.
_WAIT_POLL_SECONDS = 5.0


def _session() -> str:
    """Return this session's identity, never an anonymous one."""
    return os.environ.get(SESSION_ENV) or f"unnamed@pid{os.getpid()}"


def _build_commit(repo: Path) -> str:
    """Return HEAD, or ``"unknown"`` — never a fabricated hash."""
    try:
        done = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                              capture_output=True, text=True, timeout=15, check=False)
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    return done.stdout.strip() if done.returncode == 0 else "unknown"


def _config_hash(command: list[str]) -> str:
    """Hash the command as the run's configuration.

    The driver's own record hashes its parsed config, which is stronger; this is the pre-run
    approximation available BEFORE the driver has parsed anything. It is deliberately the whole
    argument vector minus the output path, so two runs differing only in where they write are
    correctly seen as the same measurement.

    Args:
        command: The argument vector, driver script first.

    Returns:
        ``sha256:`` hex digest.
    """
    filtered, skip = [], False
    for token in command:
        if skip:
            skip = False
            continue
        if token in ("--out", "--output", "--build-commit"):
            skip = True
            continue
        if token.startswith(("--out=", "--output=", "--build-commit=")):
            continue
        filtered.append(token)
    canonical = json.dumps(filtered, separators=(",", ":"))
    return "sha256:" + hashlib.sha256(canonical.encode()).hexdigest()


def _closure_digest(closure: list[Path], repo: Path) -> str:
    """Return a sha256 over the CONTENT of every first-party file the driver imports.

    ``build_commit`` records the repo HEAD, which changes whenever ANY file is committed — so two runs of
    byte-identical code carry different commits if an unrelated commit landed between them, and a sweep
    that compares its points by commit refuses runs that are in fact the same code. This digest answers
    the question the commit is usually being asked to answer: *did the code that ran change?*

    Path-ordered and content-only, so it is stable across machines and independent of mtime.

    Args:
        closure: Repo-relative first-party paths, as :func:`import_closure` returns them.
        repo: Repository root.

    Returns:
        ``sha256:`` hex digest, or ``""`` if any file could not be read — an incomplete digest must not
        masquerade as a complete one.
    """
    digest = hashlib.sha256()
    for path in sorted(closure):
        try:
            digest.update(str(path).encode())
            digest.update(b"\0")
            digest.update((repo / path).read_bytes())
        except OSError:
            return ""
    return "sha256:" + digest.hexdigest()


def import_closure(driver: Path, repo: Path) -> list[Path]:
    r"""Return the driver plus every first-party module it transitively imports, repo-relative.

    Resolution is by **AST**, not by regex and not by importing: the file is parsed, its ``import`` and
    ``from … import`` statements are read, and any dotted name under ``ffn_sim`` is mapped to a file
    (``a/b.py`` or ``a/b/__init__.py``) and followed.  A ``from aleph.pkg import name`` is ambiguous
    between a submodule and an attribute, so BOTH candidates are tried and whichever exists is taken.

    This deliberately does not attempt dynamic imports (``importlib``, function-local imports of
    third-party names).  Function-local imports of first-party modules ARE caught, because the walk
    visits every node rather than only module level — which matters here, since several engine modules
    import their kernels inside methods.

    Args:
        driver: The driver script, absolute or repo-relative.
        repo: Repository root.

    Returns:
        Repo-relative paths, sorted and unique, of the driver and its first-party closure.  Files that
        cannot be parsed are skipped rather than raising: the guard's job is to compare what it can
        establish, and a syntax error will surface on the run itself.
    """
    import ast

    start = driver if driver.is_absolute() else repo / driver
    pending, seen, out = [start.resolve()], set(), []
    while pending:
        path = pending.pop()
        if path in seen or not path.is_file():
            continue
        seen.add(path)
        try:
            out.append(path.relative_to(repo))
        except ValueError:                       # outside the repo — not ours to verify
            continue
        try:
            tree = ast.parse(path.read_text(encoding="utf-8", errors="replace"), filename=str(path))
        except SyntaxError:
            continue
        modules: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                modules.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
                modules.add(node.module)
                # `from pkg import name` — `name` may be a submodule or an attribute; try both.
                modules.update(f"{node.module}.{alias.name}" for alias in node.names)
        for module in modules:
            if module != "aleph" and not module.startswith("aleph."):
                continue
            base = repo.joinpath(*module.split("."))
            for candidate in (base.with_suffix(".py"), base / "__init__.py"):
                if candidate.is_file():
                    pending.append(candidate.resolve())
    return sorted(set(out))


def remote_mismatches(
    paths: list[Path], repo: Path, *, host: str, remote_root: str,
) -> tuple[list[str], list[str]]:
    r"""Compare local file contents against the remote tree the run will actually execute.

    WHY THIS EXISTS.  ``run`` stamps the LOCAL ``HEAD`` onto whatever code the remote tree happens to
    hold, and that tree is updated by hand — Syncthing shares only ``aleph/outputs/`` in this
    repository, so source has never synced by itself.  Forget one rsync and the result is a run record
    carrying a build commit that does not describe the code that ran, which ``STATE.md`` (b) names as
    the audit's top finding and which cannot be caught by reading the record afterwards.

    The comparison is on the WORKING TREE, not on ``HEAD``: the working tree is what an rsync copies
    and therefore what the remote can be made to match.  A dirty file that has been pushed is correctly
    seen as agreeing, and the run's own ``build_stamp`` records the dirty flag separately.

    Args:
        paths: Repo-relative files to compare.
        repo: Repository root.
        host: ssh host of the remote tree.
        remote_root: Root of the remote tree, as a remote-shell path.

    Returns:
        ``(mismatched, unverifiable)`` — repo-relative path strings whose remote digest differs, and
        those the remote could not be asked about (missing there, or the check itself failed).
    """
    digests: dict[str, str] = {}
    for relative in paths:
        local = repo / relative
        try:
            digests[str(relative)] = hashlib.sha256(local.read_bytes()).hexdigest()
        except OSError:
            continue
    if not digests:
        return [], []

    # Names go over stdin REPO-RELATIVE and the remote shell cds first. Sending them prefixed with
    # `remote_root` does not work when that root is `~/…`: tilde expansion is a shell feature and
    # stdin data is never expanded, so every path would arrive literal and every file would look
    # absent — which is a silent PASS-shaped failure, since "absent" and "matching" both produce no
    # mismatch unless absence is reported separately. It is, and that is how this was caught.
    listing = "\n".join(digests)
    remote_command = f"cd {remote_root} && xargs -d '\\n' -r sha256sum 2>/dev/null || true"
    try:
        done = subprocess.run(["ssh", host, remote_command], input=listing,
                              capture_output=True, text=True, timeout=120, check=False)
    except (OSError, subprocess.SubprocessError):
        return [], sorted(digests)

    remote: dict[str, str] = {}
    prefix = remote_root.rstrip("/") + "/"
    for line in done.stdout.splitlines():
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        digest, name = parts[0], parts[1].strip()
        if name.startswith(prefix):
            name = name[len(prefix):]
        remote[name] = digest

    mismatched = [name for name, digest in digests.items()
                  if name in remote and remote[name] != digest]
    unverifiable = [name for name in digests if name not in remote]
    return sorted(mismatched), sorted(unverifiable)


def _print_lease(prefix: str = "") -> int:
    lease = read_lease()
    if lease is None or not lease.is_live():
        print(f"{prefix}GPU is FREE")
        return 0
    print(f"{prefix}GPU held by {lease.holder!r} on {lease.host} (pid {lease.pid})")
    print(f"{prefix}  reason:    {lease.reason or '(unstated)'}")
    print(f"{prefix}  remaining: {lease.seconds_remaining() / 60.0:.1f} min")
    return 0


def _cmd_status(args: argparse.Namespace) -> int:
    _print_lease()
    print()
    summary = RunIndex().summary()
    print(f"run index: {summary['n_runs']} runs, {summary['n_distinct_configs']} distinct configs, "
          f"{summary['n_builds']} builds  ({summary['path']})")
    return 0


def _cmd_acquire(args: argparse.Namespace) -> int:
    try:
        lease = acquire(_session(), minutes=args.minutes, reason=args.reason)
    except LeaseHeld as held:
        print(f"REFUSED — {held}", file=sys.stderr)
        return 2
    print(f"acquired for {args.minutes:g} min by {lease.holder!r}; release when done")
    return 0


def _cmd_release(args: argparse.Namespace) -> int:
    try:
        removed = release(_session(), force=args.force)
    except LeaseHeld as held:
        print(f"REFUSED — {held}\n(the expiry will free it; --force only if the holder is known dead)",
              file=sys.stderr)
        return 2
    print("released" if removed else "nothing to release")
    return 0


def _cmd_index(args: argparse.Namespace) -> int:
    records = list(RunIndex())[-args.limit:]
    if not records:
        print("run index is empty")
        return 0
    for r in records:
        stamp = time.strftime("%m-%d %H:%M", time.localtime(r.started_at))
        print(f"{stamp}  {r.build_commit[:8]}  {r.holder:<16.16}  {r.population:<28.28}  {r.run_label}")
    return 0


def _cmd_run(args: argparse.Namespace) -> int:
    """Take the lease, refuse a duplicate, run natively, record it, and always release."""
    repo = Path(__file__).resolve().parents[2]
    command = args.command
    if not command:
        print("nothing to run", file=sys.stderr)
        return 2

    driver = repo / command[0]
    if not driver.exists():
        print(f"driver not found: {command[0]} — refusing BEFORE taking the lease so a typo does "
              f"not idle the device", file=sys.stderr)
        return 2

    build = _build_commit(repo)

    # THE BUILD STAMP HAS TO BE TRUE BEFORE THE DEVICE IS SPENT. `--build-commit` below asserts the
    # local HEAD describes the code that runs, and nothing made that so: the remote tree is updated by
    # hand. Checked here, before the lease, so a stale tree costs nothing.
    if not args.allow_stale:
        closure = import_closure(Path(command[0]), repo)
        closure_digest = _closure_digest(closure, repo)
        mismatched, unverifiable = remote_mismatches(
            closure, repo, host=args.host, remote_root=args.remote_root)
        if mismatched or unverifiable:
            print(f"REFUSED — {args.host}:{args.remote_root} does not hold the code this run would "
                  f"stamp as {build[:8]} ({len(closure)} first-party files checked):", file=sys.stderr)
            for name in mismatched[:12]:
                print(f"    differs:      {name}", file=sys.stderr)
            for name in unverifiable[:12]:
                print(f"    absent there: {name}", file=sys.stderr)
            if len(mismatched) + len(unverifiable) > 24:
                print(f"    … and {len(mismatched) + len(unverifiable) - 24} more", file=sys.stderr)
            print("  A run record whose build commit does not describe the code that ran cannot be "
                  "corrected by reading it later, which is why this refuses rather than warns.",
                  file=sys.stderr)
            print(f"  Fix:  rsync -a --relative {' '.join(str(p) for p in closure[:3])} … "
                  f"{args.host}:{args.remote_root}/", file=sys.stderr)
            print("  Or pass --allow-stale to record the run as deliberately unverified.",
                  file=sys.stderr)
            return 4
        print(f"[ffn_gpu] remote tree matches local for {len(closure)} first-party files")

    config = _config_hash(command)
    prior = RunIndex().find(build_commit=build, config_sha256=config)
    if prior and not args.force:
        print(f"REFUSED — this exact (build {build[:8]}, config {config[12:20]}) already ran "
              f"{len(prior)}x:", file=sys.stderr)
        for r in prior[-3:]:
            print(f"    {time.strftime('%m-%d %H:%M', time.localtime(r.started_at))}  "
                  f"{r.holder}  -> {r.artifact}", file=sys.stderr)
        print("    re-run anyway with --force (it is then recorded as deliberate)", file=sys.stderr)
        return 3

    reason = args.reason or " ".join(command[:1])
    try:
        acquire(_session(), minutes=args.minutes, reason=reason)
    except LeaseHeld as held:
        if not args.wait:
            print(f"REFUSED — {held}", file=sys.stderr)
            print("    --wait blocks until the device frees instead of giving up", file=sys.stderr)
            return 2
        deadline = time.time() + args.wait * 60.0
        print(f"[ffn_gpu] {held}\n[ffn_gpu] waiting up to {args.wait:g} min for the device")
        while time.time() < deadline:
            time.sleep(_WAIT_POLL_SECONDS)
            try:
                acquire(_session(), minutes=args.minutes, reason=reason)
                break
            except LeaseHeld as still:
                held = still
        else:
            print(f"GAVE UP after {args.wait:g} min — {held}", file=sys.stderr)
            return 2
        print(f"[ffn_gpu] acquired after waiting {(time.time() - (deadline - args.wait * 60.0)):.0f} s")

    # Only hand the build stamp to drivers that declare the option. Appending it unconditionally made
    # every driver that predates the flag die on `unrecognized arguments: --build-commit` AFTER the
    # lease was taken — a queued job then burns a device slot to print an argparse usage message. The
    # stamp is a convenience for the driver's own record; `RunIndex` records the build either way, so a
    # driver that cannot take it loses nothing except self-stamping, and that is recorded rather than
    # silently assumed.
    accepts_stamp = "build-commit" in (repo / command[0]).read_text(errors="replace")
    remote = (f"cd {args.remote_root} && PYTHONPATH={args.remote_root}:{args.remote_root}/ffn_sim "
              f"{args.remote_python} " + " ".join(shlex.quote(c) for c in command)
              + (f" --build-commit {build[:8]}" if accepts_stamp else ""))
    if not accepts_stamp:
        print(f"[ffn_gpu] {command[0]} declares no --build-commit; the run index still records "
              f"{build[:8]}, but this run's own artifact will not self-stamp its build")
    started = time.time()
    try:
        print(f"[ffn_gpu] {args.host}: {command[0]}  (build {build[:8]}, lease {args.minutes:g} min)")
        completed = subprocess.run(["ssh", args.host, remote], check=False)
        wrote = _artifact_exists(args.host, args.artifact) if args.artifact else None
        RunIndex().append(RunRecord(
            run_label=command[0], build_commit=build, config_sha256=config,
            artifact=args.artifact or "(not declared)", device=args.host,
            started_at=started, holder=_session(), population=args.population,
            closure_sha256=closure_digest,
        ))
        if wrote is False:
            print(f"[ffn_gpu] ⚠ the run exited {completed.returncode} but NOTHING was written to "
                  f"{args.artifact!r} on {args.host}.\n"
                  f"[ffn_gpu]   The usual cause is that --artifact was declared to this launcher while the "
                  f"driver was never told where to write (it takes its own --out).\n"
                  f"[ffn_gpu]   The index row above therefore names an artifact that does not exist; the "
                  f"device time is spent and the result is gone.", file=sys.stderr)
            return completed.returncode or 5
        return completed.returncode
    finally:
        release(_session(), force=True)
        print("[ffn_gpu] lease released")


def _artifact_exists(host: str, artifact: str) -> bool:
    """Return whether the declared artifact is actually on the run host after the run.

    WHY THIS IS CHECKED AT ALL.  ``--artifact`` tells the LAUNCHER what the run will produce; it does not
    tell the DRIVER where to write, because drivers take their own ``--out``.  Nothing connected the two,
    so a command that omitted ``--out`` burned a device slot, exited 0, and wrote an index row naming a
    file that was never created.  That happened on 2026-07-29 to a re-run whose whole purpose was to
    re-judge a failed gate — the most expensive kind of run to lose silently, because the loss looks
    exactly like a success until someone goes to read the result.

    A missing artifact is reported, not raised: the run itself may have been fine and the operator may
    simply have declared the wrong path, and turning that into an exception would lose the exit code the
    caller needs.  The exit code becomes 5 so a drainer does not record the job as cleanly done.

    Args:
        host: The run host, as passed to ssh.
        artifact: Repo-relative path the caller declared.

    Returns:
        True if present, False if absent or unverifiable — an unverifiable claim is not a verified one.
    """
    probe = f"test -s ~/ffn_ac_native/{shlex.quote(artifact)}"
    return subprocess.run(["ssh", host, probe], check=False,
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0


def main() -> int:
    """Parse and dispatch."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("status", help="who holds the GPU, and how big the run index is"
                   ).set_defaults(func=_cmd_status)

    p = sub.add_parser("acquire", help="take the lease")
    p.add_argument("--minutes", type=float, required=True,
                   help="lease duration — the length the job actually needs. Too long idles the "
                        "device if the session dies; too short expires mid-run and lets a "
                        "neighbour land on the same device")
    p.add_argument("--reason", default="", help="what the device is for, for whoever is waiting")
    p.set_defaults(func=_cmd_acquire)

    p = sub.add_parser("release", help="give the lease back")
    p.add_argument("--force", action="store_true",
                   help="release a live lease held by someone else — only when that holder is known "
                        "dead; the expiry exists to make this unnecessary")
    p.set_defaults(func=_cmd_release)

    p = sub.add_parser("index", help="what has actually run")
    p.add_argument("--limit", type=int, default=20)
    p.set_defaults(func=_cmd_index)

    p = sub.add_parser("run", help="native run with the lease enforced and the duplicate check applied")
    p.add_argument("--minutes", type=float, default=60.0)
    p.add_argument("--reason", default="")
    p.add_argument("--host", default="gbook")
    p.add_argument("--remote-root", default="~/ffn_ac_native")
    p.add_argument("--remote-python", default="~/miniconda3/envs/ffn_sim/bin/python")
    p.add_argument("--artifact", default="", help="where the run will write its record")
    p.add_argument("--population", default="",
                   help="one-line population census — recorded because the 2026-07-28 retraction "
                        # `%%`: argparse expands help through `% params`. See
                        # tests/scripts/test_argparse_help_is_formattable.
                        "turned on nobody seeing at a glance that a result came from 0.18%% of native")
    p.add_argument("--force", action="store_true", help="run even though this (build, config) ran before")
    p.add_argument("--wait", type=float, default=60.0, metavar="MIN",
                   help="if the device is busy, block up to MIN minutes and launch as soon as it frees "
                        "(default 60). Waiting is the DEFAULT because refusing was: on 2026-07-28 two "
                        "lanes wrote 'native run 0 — the device was leased all session' into their "
                        "reports while the holder sat at 5%% of device memory, and neither had any "
                        "option between 'acquire now' and 'give up'. Pass --wait 0 to get the old "
                        "fail-fast behaviour when you would rather do something else than queue")
    p.add_argument("--allow-stale", action="store_true",
                   help="launch even though the remote tree does not hold the code this run will "
                        "stamp as the build commit. Records the run as deliberately unverified; the "
                        "guard exists because a wrong build stamp cannot be detected after the fact")
    p.add_argument("command", nargs=argparse.REMAINDER,
                   help="-- then the driver path and its arguments")
    p.set_defaults(func=_cmd_run)

    args = ap.parse_args()
    if getattr(args, "command", None) and args.command and args.command[0] == "--":
        args.command = args.command[1:]
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
