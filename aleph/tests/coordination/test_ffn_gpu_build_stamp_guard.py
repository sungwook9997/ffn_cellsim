"""The guard that makes a run's build commit true before the device is spent.

The launcher passes the LOCAL ``HEAD`` to the driver as ``--build-commit``, and the driver stamps
it into the run record.  Nothing made that assertion true: Syncthing shares only ``aleph/outputs/`` in
this repository, so the remote source tree is updated by hand and one forgotten rsync produces a record
whose build commit does not describe the code that ran.  ``STATE.md`` (b) names a missing/wrong build
stamp as the whole-repo audit's top finding, and it is not detectable by reading the record afterwards —
which is why the check has to happen before launch, and why it refuses rather than warns.

Measured while writing these: the first version sent remote paths prefixed with ``~/ffn_ac_native``.
Tilde expansion is a shell feature and stdin data is never expanded, so every file arrived literal and
every file looked ABSENT.  That is a PASS-shaped failure unless absence is reported separately from
mismatch — it is, and that is how it was caught.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from aleph.scripts.run_provenance import import_closure, remote_mismatches


def _write(root: Path, relative: str, text: str) -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture()
def repo(tmp_path: Path) -> Path:
    """A miniature repository with the same package shape the real closure walks."""
    _write(tmp_path, "aleph/__init__.py", "")
    _write(tmp_path, "aleph/ac/__init__.py", "")
    _write(tmp_path, "aleph/ac/leaf.py", "VALUE = 1\n")
    _write(tmp_path, "aleph/ac/mid.py",
           "import numpy as np\nfrom aleph.ac.leaf import VALUE\n")
    _write(tmp_path, "aleph/ac/deep/__init__.py", "from aleph.ac.leaf import VALUE\n")
    _write(tmp_path, "aleph/scripts/driver.py",
           "import warp as wp\n"
           "from aleph.ac.mid import VALUE\n"
           "def run():\n"
           "    from aleph.ac import deep\n"          # function-local, and a package
           "    return deep\n")
    return tmp_path


def test_the_closure_follows_first_party_imports_transitively(repo: Path) -> None:
    closure = import_closure(Path("aleph/scripts/driver.py"), repo)
    names = {str(path) for path in closure}
    assert "aleph/scripts/driver.py" in names
    assert "aleph/ac/mid.py" in names                 # direct
    assert "aleph/ac/leaf.py" in names                # transitive, through mid
    assert "aleph/ac/deep/__init__.py" in names       # a PACKAGE, resolved to its __init__


def test_a_function_local_import_is_followed_because_engine_modules_use_them(repo: Path) -> None:
    """``ac/cell/implicit_mechanics`` imports its projector kernel inside a method; module-level-only
    walking would silently drop exactly the files a physics run depends on."""
    closure = import_closure(Path("aleph/scripts/driver.py"), repo)
    assert Path("aleph/ac/deep/__init__.py") in closure


def test_third_party_imports_are_not_claimed_as_ours(repo: Path) -> None:
    """`numpy` and `warp` live in the environment, not the tree; asserting on them would be noise."""
    closure = import_closure(Path("aleph/scripts/driver.py"), repo)
    assert not any("numpy" in str(path) or "warp" in str(path) for path in closure)


def test_an_unparseable_module_is_skipped_rather_than_raising(repo: Path) -> None:
    """The guard compares what it can establish; a syntax error surfaces on the run itself."""
    _write(repo, "aleph/ac/broken.py", "def (:\n")
    _write(repo, "aleph/scripts/driver2.py", "from aleph.ac.broken import thing\n")
    closure = import_closure(Path("aleph/scripts/driver2.py"), repo)
    assert Path("aleph/ac/broken.py") in closure      # still compared, just not walked into


def _fake_ssh(monkeypatch: pytest.MonkeyPatch, stdout: str) -> dict[str, str]:
    """Replace the ssh call and capture what would have been sent."""
    captured: dict[str, str] = {}

    def fake_run(argv, **kwargs):  # noqa: ANN001, ANN003
        captured["argv"] = " ".join(argv)
        captured["stdin"] = kwargs.get("input", "")
        return subprocess.CompletedProcess(argv, 0, stdout=stdout, stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    return captured


def test_matching_content_reports_no_mismatch_and_no_absence(
    repo: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    import hashlib

    relative = Path("aleph/ac/leaf.py")
    digest = hashlib.sha256((repo / relative).read_bytes()).hexdigest()
    captured = _fake_ssh(monkeypatch, f"{digest}  {relative}\n")

    mismatched, unverifiable = remote_mismatches(
        [relative], repo, host="somehost", remote_root="~/tree")
    assert mismatched == [] and unverifiable == []
    # REGRESSION: names must go over stdin repo-relative, because `~` is never expanded in stdin data
    assert captured["stdin"] == str(relative)
    assert "~" not in captured["stdin"]
    assert "cd ~/tree" in captured["argv"]


def test_different_content_is_reported_as_a_mismatch(
    repo: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    relative = Path("aleph/ac/leaf.py")
    _fake_ssh(monkeypatch, f"{'0' * 64}  {relative}\n")
    mismatched, unverifiable = remote_mismatches(
        [relative], repo, host="somehost", remote_root="~/tree")
    assert mismatched == [str(relative)]
    assert unverifiable == []


def test_a_file_the_remote_never_answered_for_is_ABSENT_not_matching(
    repo: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The failure mode that hid the tilde bug: silence must never read as agreement."""
    relative = Path("aleph/ac/leaf.py")
    _fake_ssh(monkeypatch, "")                    # remote said nothing about anything
    mismatched, unverifiable = remote_mismatches(
        [relative], repo, host="somehost", remote_root="~/tree")
    assert mismatched == []
    assert unverifiable == [str(relative)]


def test_an_absolute_prefix_in_the_remote_reply_is_stripped_before_comparing(
    repo: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    import hashlib

    relative = Path("aleph/ac/leaf.py")
    digest = hashlib.sha256((repo / relative).read_bytes()).hexdigest()
    _fake_ssh(monkeypatch, f"{digest}  /home/someone/tree/{relative}\n")
    mismatched, unverifiable = remote_mismatches(
        [relative], repo, host="somehost", remote_root="/home/someone/tree")
    assert mismatched == [] and unverifiable == []


def test_an_ssh_failure_is_unverifiable_never_silently_verified(
    repo: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    def blow_up(argv, **kwargs):  # noqa: ANN001, ANN003
        raise OSError("ssh: connection refused")

    monkeypatch.setattr(subprocess, "run", blow_up)
    mismatched, unverifiable = remote_mismatches(
        [Path("aleph/ac/leaf.py")], repo, host="somehost", remote_root="~/tree")
    assert mismatched == []
    assert unverifiable == ["aleph/ac/leaf.py"]
