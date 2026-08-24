"""What build a run was measured on, and a canonical hash for the configuration it ran.

Two jobs, both about making a record answer "which code produced this?" without the reader having to
trust anybody.

**The build stamp.**  The audit that motivated Project Aleph found that every GPU run record in the
reference project stamps its build source as ``declared`` — the GPU host is a synced tree, not a git
checkout, so the commit is the caller's word.  That is not a defect in itself; a declared commit
honestly labelled is far better than a fabricated one.  The defect is what happens when git *is*
also available and disagrees.  The reference keeps git's answer in the ``commit`` field and files the
disagreement in a side key, so a reader who looks at ``commit`` never learns it was contested.

Aleph refuses to pick.  :class:`BuildStamp` keeps ``git_commit`` and ``declared_commit`` as separate
fields, and the :attr:`BuildStamp.commit` property **raises** when they disagree.  Code that wants a
single commit must decide what to do about the contest in the open.  Serialization writes
``commit: null`` plus a ``disagreement`` block carrying both values.

**The config hash.**  :func:`config_hash` is deliberately not ``json.dumps(..., default=str)``.  A
``default=str`` fallback turns every object the encoder does not understand into its ``repr``, which
means two genuinely different configs can serialize identically and two identical configs can
serialize differently across Python versions.  Aleph's canonicalization refuses what it does not
understand.  The full rule set is documented on :func:`canonical_json`.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import unicodedata
from dataclasses import dataclass
from enum import Enum, unique
from pathlib import Path
from typing import Any, Mapping, Sequence

__all__ = [
    "BuildSource",
    "BuildStamp",
    "CommitDisagreementError",
    "RepositoryState",
    "build_stamp",
    "canonical_json",
    "config_hash",
]

#: How long a git subprocess may run before the stamp gives up and records the build as unavailable.
#: A hung git call must never hang a run; an unavailable build is a recordable fact, a hang is not.
_GIT_TIMEOUT_SECONDS: float = 15.0


class CommitDisagreementError(RuntimeError):
    """Raised when a single commit is demanded from a stamp whose two sources disagree.

    Deliberately not a warning.  The whole point of recording a disagreement instead of resolving it
    is that downstream code must not be able to proceed as though one answer had won.
    """


@unique
class BuildSource(Enum):
    """How the stamp learned which commit the run was made on.

    Members:
        GIT: A git checkout was interrogated and answered.  The commit is verified against the tree
            the stamp was taken from.
        DECLARED: The caller asserted the commit because the executing host is not a git checkout.
            Aleph did not verify it.  This is honest and usable; it is not proof.
        UNAVAILABLE: No commit could be established at all — git absent, the directory is not a
            repository, or the repository has no commits yet — and none was declared.
    """

    GIT = "git"
    DECLARED = "declared"
    UNAVAILABLE = "unavailable"


@unique
class RepositoryState(Enum):
    """What the stamp found where it looked for a repository.

    Kept separate from :class:`BuildSource` because the two answer different questions.  ``source``
    says where the commit string came from; ``repository_state`` says what was actually on disk, and
    a reader needs both to understand an ``UNAVAILABLE`` stamp.

    Members:
        NOT_A_REPOSITORY: The path is not inside a git working tree, or git is not installed.
        UNBORN: A repository exists but has no commits yet — ``HEAD`` points at a branch that does
            not exist.  This is the state of a freshly initialised repo, and it must be recorded as
            its own fact rather than collapsed into "not a repository", because the two call for
            entirely different fixes.
        COMMITTED: A repository exists and ``HEAD`` resolves to a commit.
    """

    NOT_A_REPOSITORY = "not_a_repository"
    UNBORN = "unborn"
    COMMITTED = "committed"


@dataclass(frozen=True, slots=True)
class BuildStamp:
    """The build a run was measured on, with the two possible commit sources kept apart.

    Attributes:
        source: How the commit was obtained.  See :class:`BuildSource`.
        repository_state: What was found on disk where the stamp looked.
        git_commit: The commit git reported, or ``None`` if git could not report one.
        declared_commit: The commit the caller asserted, or ``None`` if none was asserted.
        dirty: Whether the working tree had uncommitted changes.  ``None`` when unknowable.  A dirty
            tree is recorded, not refused — development runs are legitimate — but it means the
            commit alone does not reproduce the run.
        branch: The branch name, or ``None``.  Present even for an unborn repository, where it names
            the branch the first commit will land on.
        reason: Why the stamp is weaker than a verified commit, when it is.
    """

    source: BuildSource
    repository_state: RepositoryState
    git_commit: str | None = None
    declared_commit: str | None = None
    dirty: bool | None = None
    branch: str | None = None
    reason: str | None = None

    @property
    def has_disagreement(self) -> bool:
        """Whether git and the caller named different commits.

        A declared commit is treated as an abbreviation of the git commit when the git commit starts
        with it, since callers routinely declare short hashes.  Anything else is a disagreement.
        """
        if self.git_commit is None or self.declared_commit is None:
            return False
        return not self.git_commit.startswith(self.declared_commit)

    @property
    def commit(self) -> str | None:
        """The single commit this stamp attests to.

        Returns:
            The git commit when git answered, the declared commit when it did not, and ``None``
            when neither is available.

        Raises:
            CommitDisagreementError: When git and the caller disagree.  Aleph does not resolve the
                contest, so there is no single commit to return and callers must handle both values
                explicitly via :attr:`git_commit` and :attr:`declared_commit`.
        """
        if self.has_disagreement:
            raise CommitDisagreementError(
                f"this stamp has two contradictory commits — git says {self.git_commit!r}, the "
                f"caller declared {self.declared_commit!r}. Aleph records the disagreement rather "
                "than resolving it; read git_commit and declared_commit and decide in the open"
            )
        if self.git_commit is not None:
            return self.git_commit
        return self.declared_commit

    @property
    def is_reproducible(self) -> bool:
        """Whether this stamp alone identifies the code that ran.

        ``False`` for a dirty tree, a contested commit, a declared-only commit, or no commit at all.
        A record carrying ``is_reproducible: False`` is still a valid record; it just cannot be the
        sole basis for a reproduction claim.
        """
        return (
            self.source is BuildSource.GIT
            and not self.has_disagreement
            and self.dirty is False
            and self.git_commit is not None
        )

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able build block.

        A contested stamp writes ``commit: null`` and a ``disagreement`` block.  Both raw values stay
        in their own fields regardless, so no information is lost by the contest.
        """
        block: dict[str, Any] = {
            "source": self.source.value,
            "repository_state": self.repository_state.value,
            "commit": None if self.has_disagreement else self.commit,
            "git_commit": self.git_commit,
            "declared_commit": self.declared_commit,
            "dirty": self.dirty,
            "branch": self.branch,
            "is_reproducible": self.is_reproducible,
        }
        if self.reason:
            block["reason"] = self.reason
        if self.has_disagreement:
            block["disagreement"] = {
                "git": self.git_commit,
                "declared": self.declared_commit,
                "resolved": False,
                "note": (
                    "the executing host declared a commit that the git tree does not confirm. "
                    "Aleph records both and resolves neither: preferring either one would present "
                    "a contested build as a settled one"
                ),
            }
        return block


def _run_git(root: Path, *args: str) -> str | None:
    """Run one git command in ``root`` and return trimmed stdout, or ``None`` on any failure.

    Args:
        root: Directory to run git in.
        *args: Git arguments after ``-C <root>``.

    Returns:
        Stripped stdout on exit status zero; ``None`` on non-zero status, a missing git binary, or a
        timeout.  Failures are collapsed to ``None`` on purpose: the caller distinguishes the cases
        it cares about by *which* commands failed, not by their error text.
    """
    try:
        completed = subprocess.run(
            ("git", "-C", str(root), *args),
            capture_output=True,
            text=True,
            timeout=_GIT_TIMEOUT_SECONDS,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode != 0:
        return None
    return completed.stdout.strip()


def build_stamp(
    repo_root: Path | str | None = None,
    *,
    declared_commit: str | None = None,
) -> BuildStamp:
    """Capture the build a run is being measured on.

    Four outcomes, each recorded as a distinct state rather than blurred together:

    * **git answered.** ``source=GIT``, ``repository_state=COMMITTED``.  If a commit was also
      declared and it disagrees, both are kept and :attr:`BuildStamp.commit` will raise.
    * **the repository exists but has no commits.** ``repository_state=UNBORN``.  A freshly
      initialised repo is not the same thing as no repo, and reporting it as "not a repository"
      would send someone looking for a missing checkout instead of a missing commit.  The branch and
      dirty flag are still recorded, because both are knowable and both are useful.
    * **no repository, but a commit was declared.** ``source=DECLARED``.  This is the normal state
      for a run on a synced tree or a container image.
    * **nothing at all.** ``source=UNAVAILABLE``, with a reason.  Never a fabricated hash, never a
      placeholder that reads like a commit.

    Args:
        repo_root: Directory to resolve the repository from.  Defaults to the Aleph repository this
            module lives in.
        declared_commit: The commit the caller asserts the run was made on, for hosts that are not
            git checkouts.  Recorded as a declaration and cross-checked when git is also available.

    Returns:
        The stamp.

    Raises:
        ValueError: If ``declared_commit`` is supplied but is not a plausible hex object name.  A
            declared commit that cannot be a commit is a typo, and recording it would put a
            meaningless string where a build identity belongs.
    """
    if declared_commit is not None:
        cleaned = declared_commit.strip().lower()
        if len(cleaned) < 7 or len(cleaned) > 64 or any(c not in "0123456789abcdef" for c in cleaned):
            raise ValueError(
                f"declared_commit must be 7-64 hex characters; got {declared_commit!r}"
            )
        declared_commit = cleaned

    root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[2]

    inside_tree = _run_git(root, "rev-parse", "--is-inside-work-tree") == "true"
    head_commit = _run_git(root, "rev-parse", "--verify", "HEAD") if inside_tree else None

    if inside_tree and head_commit:
        porcelain = _run_git(root, "status", "--porcelain")
        return BuildStamp(
            source=BuildSource.GIT,
            repository_state=RepositoryState.COMMITTED,
            git_commit=head_commit,
            declared_commit=declared_commit,
            dirty=None if porcelain is None else bool(porcelain),
            branch=_run_git(root, "rev-parse", "--abbrev-ref", "HEAD"),
        )

    if inside_tree:
        # A repository with an unborn HEAD. `rev-parse --abbrev-ref HEAD` fails here, but
        # `symbolic-ref` still names the branch the first commit will create.
        porcelain = _run_git(root, "status", "--porcelain")
        branch = _run_git(root, "symbolic-ref", "--short", "HEAD")
        unborn_reason = (
            f"{root} is a git repository with no commits yet, so there is no object name to stamp. "
            "Nothing here is reproducible until the first commit lands"
        )
        if declared_commit is not None:
            return BuildStamp(
                source=BuildSource.DECLARED,
                repository_state=RepositoryState.UNBORN,
                git_commit=None,
                declared_commit=declared_commit,
                dirty=None if porcelain is None else bool(porcelain),
                branch=branch,
                reason=(
                    unborn_reason
                    + "; the recorded commit is the caller's assertion and was not verified here"
                ),
            )
        return BuildStamp(
            source=BuildSource.UNAVAILABLE,
            repository_state=RepositoryState.UNBORN,
            dirty=None if porcelain is None else bool(porcelain),
            branch=branch,
            reason=unborn_reason,
        )

    if declared_commit is not None:
        return BuildStamp(
            source=BuildSource.DECLARED,
            repository_state=RepositoryState.NOT_A_REPOSITORY,
            declared_commit=declared_commit,
            reason=(
                f"{root} is not a git working tree, so the commit is the caller's assertion and "
                "was not verified here"
            ),
        )

    return BuildStamp(
        source=BuildSource.UNAVAILABLE,
        repository_state=RepositoryState.NOT_A_REPOSITORY,
        reason=(
            f"git is unavailable or {root} is not a working tree, and no commit was declared; the "
            "build this run was measured on cannot be established"
        ),
    )


def _canonical_key(key: object) -> str:
    """Return the canonical text form of a mapping key.

    Args:
        key: A mapping key.

    Returns:
        The NFC-normalized key.

    Raises:
        TypeError: If the key is not a string.  Integer and tuple keys are refused rather than
            stringified, because ``{1: "a"}`` and ``{"1": "a"}`` are different configurations and a
            hash that cannot tell them apart is not a config hash.
    """
    if not isinstance(key, str):
        raise TypeError(
            f"config keys must be strings for canonical ordering; got {type(key).__name__} "
            f"({key!r}). Convert it deliberately at the call site so the choice is visible"
        )
    return unicodedata.normalize("NFC", key)


def _encode(value: object) -> str:
    """Encode one value into Aleph's canonical JSON text.  See :func:`canonical_json` for the rules.

    Args:
        value: The value to encode.

    Returns:
        A JSON fragment.

    Raises:
        TypeError: On a type with no canonical form.
        ValueError: On a non-finite float, or on duplicate keys after NFC normalization.
    """
    if value is None:
        return "null"
    # bool before int: bool is a subclass of int, and True must not encode as 1.
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if value != value or value in (float("inf"), float("-inf")):
            raise ValueError(
                f"config values must be finite; got {value!r}. NaN and infinity have no portable "
                "JSON spelling, and a config that contains one is not reproducible from its record"
            )
        # repr gives the shortest string that round-trips to the same double, and always carries a
        # '.' or an exponent, so a float never collides with the int of the same value. float() is
        # applied first because numpy's float64 subclasses float but reprs as `np.float64(1.0)`.
        return repr(float(value))
    if isinstance(value, str):
        return json.dumps(unicodedata.normalize("NFC", value), ensure_ascii=False)
    if isinstance(value, (bytes, bytearray)):
        return json.dumps("base16:" + bytes(value).hex())
    if isinstance(value, Mapping):
        items = []
        seen: set[str] = set()
        for raw_key, item in value.items():
            key = _canonical_key(raw_key)
            if key in seen:
                raise ValueError(
                    f"key {key!r} appears twice after Unicode NFC normalization; two visually "
                    "distinct keys collapsing to one makes the hash depend on iteration order"
                )
            seen.add(key)
            items.append((key, item))
        items.sort(key=lambda pair: pair[0].encode("utf-8"))
        body = ",".join(f"{json.dumps(k, ensure_ascii=False)}:{_encode(v)}" for k, v in items)
        return "{" + body + "}"
    if isinstance(value, (set, frozenset)):
        raise TypeError(
            "sets have no canonical order that means anything to a reader; convert to a sorted "
            "list at the call site so the record shows the order that was intended"
        )
    if isinstance(value, Sequence):
        return "[" + ",".join(_encode(item) for item in value) + "]"
    # numpy scalars and arrays expose tolist(); using it keeps numpy out of this module's imports
    # while still canonicalizing the values a config realistically carries.
    tolist = getattr(value, "tolist", None)
    if callable(tolist):
        return _encode(tolist())
    raise TypeError(
        f"no canonical form for {type(value).__name__} ({value!r}). Aleph refuses to fall back to "
        "repr(): a repr fallback lets two different configs hash the same and lets one config hash "
        "differently across Python versions"
    )


def canonical_json(obj: object) -> str:
    """Return Aleph's canonical JSON text for ``obj``.

    The canonicalization, in full, so a reader can reproduce a hash by hand:

    * **Key ordering.** Mapping keys must be strings.  They are Unicode NFC-normalized, then sorted
      by their UTF-8 byte sequence.  Two keys that normalize to the same string are an error, not a
      silent overwrite.
    * **Non-string keys.** Refused.  ``{1: "a"}`` and ``{"1": "a"}`` are different configs.
    * **None.** Encodes as ``null``.  A key present with value ``None`` and a key absent entirely
      produce *different* text, and therefore different hashes.  Aleph does not strip nulls: "this
      knob was explicitly unset" and "this knob does not exist" are different configurations.
    * **Booleans.** ``true`` / ``false``, checked before integers so ``True`` never encodes as ``1``.
    * **Integers.** Plain decimal, arbitrary precision, no exponent.
    * **Floats.** ``repr()``, which is the shortest decimal that round-trips to the identical
      double.  The result always contains ``.`` or ``e``, so ``1`` and ``1.0`` never collide.
      ``-0.0`` encodes as ``-0.0`` and is therefore distinct from ``0.0``, because a sign on zero
      can be a meaningful direction convention.  NaN and infinity are refused.
    * **Strings.** NFC-normalized, then JSON-escaped with ``ensure_ascii=False`` so the text is
      UTF-8 and stable.
    * **Bytes.** Encoded as the string ``"base16:<hex>"``.
    * **Sequences.** Order-preserving, since order in a list is information.
    * **Sets.** Refused; sort at the call site so the record shows an intended order.
    * **Anything else.** Refused if it has no ``tolist()``.  There is no ``repr`` fallback.
    * **Separators.** ``,`` and ``:`` with no whitespace anywhere.

    Args:
        obj: The object to canonicalize, usually a config mapping.

    Returns:
        The canonical text.  Parsing it with :func:`json.loads` yields a structure equal to the
        input for every type JSON can represent.

    Raises:
        TypeError: On a type with no canonical form, or a non-string mapping key.
        ValueError: On a non-finite float or a post-normalization duplicate key.
    """
    return _encode(obj)


def config_hash(obj: object) -> str:
    """Return the SHA-256 of ``obj``'s canonical JSON, prefixed with its algorithm.

    Args:
        obj: The object to hash, usually a config mapping.

    Returns:
        ``"sha256:"`` followed by 64 lowercase hex characters.  The prefix is part of the value so
        that a stored digest carries its own algorithm and a future migration to a different hash
        is visible in every record rather than inferred from the length.

    Raises:
        TypeError: Propagated from :func:`canonical_json`.
        ValueError: Propagated from :func:`canonical_json`.
    """
    text = canonical_json(obj)
    return "sha256:" + hashlib.sha256(text.encode("utf-8")).hexdigest()
