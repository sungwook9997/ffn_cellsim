#!/usr/bin/env python
r"""Run provenance: does the code on the run host match the commit the record will stamp?

Split out of `ffn_gpu.py` by decision **B2**, which said to delete that script. It could not be
deleted whole, and the reason is the useful part of this module's existence.

**What was deleted**: the lease. One session took a lock on one A5000 before launching, and the
lock lived on this machine because every contender did. All three premises are gone — the A5000
was retired 2026-08-04, the device is now a Slurm-scheduled workstation shared with other people,
and `gpu-submit` refuses a busy card. Outside a Slurm allocation `cuInit()` returns
`CUDA_ERROR_NO_DEVICE`, so the operating system holds the door the lease was bolted to. Worse than
redundant: it expired on wall-clock rather than on work, which cut unattended runs, and it defaulted
to a host that no longer exists.

**What survived, and matters more now than it did then**: a run record stamps a build commit, and
`remote_mismatches` is what checks that the tree on the run host actually holds that code. A record
whose commit does not describe what ran cannot be corrected afterwards — the run is over and the
number is in a file. With the device now shared with people outside this project, "the host had what
I thought it had" is exactly the assumption that stops being safe.

`import_closure` is how that check knows what to compare: the driver's first-party import closure,
not the whole tree, so a stale file nobody imports does not fail a run and a stale file the driver
does import cannot pass one.

Nothing here takes a lock, launches anything, or reaches a device. Submit GPU work with
`MEM_GB=<GB> CPUS_PER_TASK=<n> gpu-submit <card> <time> <command>`, and ask the PI first —
`CLAUDE.md` §Running on the GPU.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

#: How a session names itself. Falls back to something a human can still trace.
SESSION_ENV = "FFN_SESSION"


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
        ``sha256:`` hex digest, or ``""`` if the closure is EMPTY or any file could not be read — an
        incomplete digest must not masquerade as a complete one, and a digest over nothing is the
        emptiest form of that.
    """
    if not closure:
        # ⚠ An EMPTY closure must not produce a digest. sha256 of nothing is a perfectly well-formed
        # hex string, so every driver whose closure failed to resolve would stamp the SAME value and
        # compare equal to every other — "this is the same code" asserted from having read none of
        # it. Found 2026-08-21 by the Sanity Gate below, in code that had been in the tree for weeks.
        return ""
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


def _runtime() -> dict[str, str]:
    """Return what compiled and ran the kernels, never raising.

    Separate from :func:`stamp`'s body so a driver in an environment without Warp still gets a
    provenance block rather than an exception — the point of the block is to say what it can, and
    ``"absent"`` is a true answer that a missing key is not.
    """
    out = {"python": sys.version.split()[0], "env_prefix": sys.prefix}
    try:
        import warp                                   # noqa: PLC0415 — probe, not a dependency
        out["warp"] = getattr(warp, "__version__", "unknown")
    except Exception:                                 # noqa: BLE001 — see docstring
        out["warp"] = "absent"
    return out


def stamp(driver: Path | str, repo: Path | None = None) -> dict[str, str]:
    """Return the provenance block a run record must carry, computed WHERE THE RUN RUNS.

    **Why this is a separate entry point from :func:`remote_mismatches`.** That function is the
    pre-flight check, and it can only run somewhere that has both trees. This one runs on the run
    host, from inside the driver, and answers the question the record has to answer by itself:
    *which code produced this number?*

    On 2026-08-21 that gap was measured rather than argued. The GPU host's tree is deployed BY HAND
    — the module docstring above and :func:`remote_mismatches` both say so — and a night's worth of
    fixes to ``step.py``, ``geometry.py`` and ``bond.py`` were never pushed. Records written after
    that point named no commit at all, so nothing in the artifact could show it; the drift was found
    by hashing the remote tree three hours later, and could just as easily not have been. **The guard
    that would have caught it existed, was tested, and was called by no production driver.**

    ``closure_digest`` is the load-bearing field, not ``build_commit``: the host has no ``.git``, so
    the commit is read from a ``COMMIT_STAMP`` file that is only as fresh as the last hand-deploy —
    exactly the thing that went stale. The digest is computed from the bytes that were imported, so
    it cannot go stale, and comparing it against the same digest taken at a local commit is what
    makes a record attributable after the fact.

    Args:
        driver: The driver script whose import closure defines "the code that ran".
        repo: Repository root. Defaults to the root inferred from this file's location.

    Returns:
        Mapping with ``build_commit`` (or ``"unknown"``), ``closure_digest`` (``""`` if incomplete),
        ``closure_files`` and ``session``. Never raises: a driver must not lose a completed run to
        its own bookkeeping.
    """
    root = Path(repo) if repo is not None else Path(__file__).resolve().parents[2]
    try:
        closure = import_closure(Path(driver), root)
        commit = _build_commit(root)
        if commit == "unknown":                       # no .git — the deployed tree. Fall back to the
            stamp_file = root / "COMMIT_STAMP"        # hand-written stamp, LABELLED as hand-written.
            if stamp_file.is_file():
                commit = f"{stamp_file.read_text().strip()} (COMMIT_STAMP, hand-deployed)"
        return {"build_commit": commit, "closure_digest": _closure_digest(closure, root),
                "closure_files": str(len(closure)), "session": _session(),
                # ⚠ The invocation, because a record that does not carry it is only reproducible by
                # someone who still has the shell it was typed in. Found the same morning as the
                # stale stamp, reconstructing a run's --bind flags out of its own output.
                "argv": " ".join(sys.argv),
                # ⚠ And the RUNTIME, because the same morning found TWO Warp versions on the run host
                # — production `ffn_sim` at 1.14.0 and `aleph` at 1.15.0 — with nothing in any record
                # saying which produced it. `closure_digest` fixes which source ran; this fixes what
                # compiled it, and for a device kernel the compiler is not a detail.
                **_runtime()}
    except Exception as exc:                          # noqa: BLE001 — see docstring: never lose a run
        return {"build_commit": "unknown", "closure_digest": "",
                "closure_files": "0", "session": f"stamp failed: {exc!r}"}


def _demo() -> None:
    """Sanity Gate for :func:`stamp` — the properties a record's provenance has to have.

    What is checkable without a device: that the digest tracks CONTENT and not the clock, that an
    unreadable closure yields ``""`` rather than a short digest that would compare equal to nothing,
    and that the function cannot take a run down with it.
    """
    import tempfile

    root = Path(__file__).resolve().parents[2]
    here = Path("aleph/scripts/run_provenance.py")

    # 1. It returns the four fields a record needs, and the digest is non-empty for a real closure.
    s = stamp(here, root)
    assert set(s) == {"build_commit", "closure_digest", "closure_files", "session", "argv",
                      "python", "env_prefix", "warp"}, s
    # ⚠ The runtime probe must answer even where Warp is absent — "absent" is a true answer and a
    # missing key is not. The dev machine has no CUDA but does have Warp, so this asserts the shape
    # rather than the value.
    assert s["warp"] and s["python"] and s["env_prefix"]
    assert s["closure_digest"].startswith("sha256:"), s["closure_digest"]
    assert int(s["closure_files"]) >= 1

    # 2. Stable across calls — content-only, so nothing about the clock or mtime leaks in.
    assert stamp(here, root)["closure_digest"] == s["closure_digest"]

    # 3. The failure mode that matters: an INCOMPLETE digest must not look like a complete one.
    #    A closure naming a file that is not there yields "", never a digest over what happened to
    #    be readable — that is the difference between "cannot attribute" and "attributed wrongly".
    assert _closure_digest([Path("aleph/definitely_absent_xyz.py")], root) == ""

    # 4. A driver must never lose a finished run to its own bookkeeping.
    assert stamp(Path(tempfile.gettempdir()) / "not_a_driver_xyz.py")["closure_digest"] == ""

    # 5. Content, not path: two different files cannot collide, and the same bytes under a different
    #    name must NOT compare equal — the path is hashed in for exactly that reason.
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        (tmp / "a.py").write_text("x = 1\n")
        (tmp / "b.py").write_text("x = 1\n")
        assert _closure_digest([Path("a.py")], tmp) != _closure_digest([Path("b.py")], tmp)
        (tmp / "a.py").write_text("x = 2\n")
        assert _closure_digest([Path("a.py")], tmp) != _closure_digest([Path("b.py")], tmp)

    print(f"run_provenance self-check OK — {s['closure_files']} files, {s['closure_digest'][:23]}…")


if __name__ == "__main__":
    _demo()
