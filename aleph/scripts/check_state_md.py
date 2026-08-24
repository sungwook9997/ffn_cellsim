"""Structural gate for the root ``STATE.md``.

``STATE.md`` is the one status document the charter permits editing at every milestone. Its value
depends entirely on three properties that a human reviewer will not re-check by hand:

1. **Every tier-(a) result row carries a build commit that exists in git.** This is the direct guard
   against the 2026-07-25 whole-repo audit's highest-severity finding — the gamma-floor headline was
   measured on a cortex build the project later declared broken and *nothing recorded which build*.
   A row with a missing, placeholder, or fabricated build commit fails here.
2. **Every path in the minimal read set exists on disk.** A read set that points at a deleted or
   renamed document silently stops being a read set.
3. **The file stays short and keeps its six sections.** Length is a declared proxy for the "~2 pages,
   nothing else" rule; the section check stops a rewrite from dropping the non-quotable list.

What this gate deliberately does NOT do (do not read it as more than it is):

- It cannot tell whether a listed result is *true*, whether its caveats are complete, or whether a
  build commit is the *right* one for that measurement. It checks that a build commit is named and
  resolves — not that the run was made from it. No artifact in the repo stamps its own build commit;
  until one does, that binding is human.
- It does not check the (c) non-quotable list for completeness, and it cannot notice a result that was
  quietly moved out of (c) into (b).
- The SIZE assertion is a proxy. A budget of confidently wrong status would pass it.
- It used to count LINES, and that was worse than a proxy — it was defeatable without deleting a word.
  On 2026-07-28 a session brought this file under the cap by re-wrapping the prose at a wider column;
  nothing was removed and the check went green. STATE.md was 225 bytes per line against CLAUDE.md's 88,
  so "140 lines" was asserting a page count roughly 2.6x wrong in the direction that hurt. The budget is
  now BYTES, which is what a reader and a context window both actually pay.

Usage:
    python aleph/scripts/check_state_md.py [--file PATH]

Exit 0 = PASS, exit 1 = FAIL with one line per finding.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_STATE = REPO_ROOT / "STATE.md"

#: Size budget in BYTES, not lines — see the module docstring for why lines were unenforceable.
#: 20,000 B is ~5,000 tokens. The file was 32,454 B on 2026-07-28 morning; moving (b)'s evidence to
#: `docs/state_rows.yaml`, (d)'s reasoning to `docs/READ_SET.md` and (f)'s narrative to
#: `docs/STATE_LOG.md` took it to ~18,000, and the budget is set just above that rather than at
#: whatever the file happens to weigh — a budget set to today's size is loosening by default.
#: Lowering it is the intended direction; raising it is a PI decision, not a convenience.
#:
#: RAISED 20,000 -> 22,000 by PI decision, 2026-08-10 (session `47e953f4-93d6-4c67-a901-08ca181222f6`).
#: Why it was asked rather than taken: the file stood at 19,959 B — 41 B of headroom — and a sixteenth
#: tier-(a) row (the Card-5 native array-ownership gate) does not fit at any honest length. Editing a
#: threshold so a gate passes is what the charter forbids, so the session halted and surfaced instead of
#: trimming to fit. The 2,000 B is deliberately small: it buys roughly four rows, not a new narrative
#: budget, and the next session to hit this ceiling should reach for (b)'s evidence layer again before
#: reaching for this number. This is the second budget to hit its ceiling in one night — see (e) 3, the
#: open boot-budget item — which is the signal worth reading, not the constant.
#: SCOPED to everything EXCEPT section (c) by PI decision, 2026-08-21; the number was NOT raised.
#: The two rules had collided: (c)'s own header says it "only GROWS; never shorten it by dropping an
#: entry", and a monotonically growing list cannot live under a fixed cap -- only the DATE of the
#: collision is negotiable. It arrived 2026-08-21: 39 B free, next row 141 B, and `c19` sat in the
#: companion file UNPUSHED, which is the one thing (c) exists to prevent. Raising the cap a second
#: time defers the same collision; exempting (c) removes it, and leaves the cap doing what it was
#: for -- stopping (b) and the prose from growing. (b) alone was 11,597 B of the 21,961.
#: The exemption is ANNOUNCED on every run (both numbers printed): an invisible carve-out stops
#: being a proxy and becomes a hole.
MAX_BYTES = 22_000

#: Section (c) is exempt from MAX_BYTES; these two headings fence it. REQUIRED_SECTIONS already
#: asserts each appears exactly once, so no separate fence markers are needed.
BLOCKLIST_BEGIN = "## (c)"
BLOCKLIST_END = "## (d)"
REQUIRED_SECTIONS = ("## (a)", "## (b)", "## (c)", "## (d)", "## (e)", "## (f)")
TIER_A_BEGIN = "<!-- STATE-TIER-A:BEGIN -->"
TIER_A_END = "<!-- STATE-TIER-A:END -->"
READSET_BEGIN = "<!-- STATE-READSET:BEGIN -->"
READSET_END = "<!-- STATE-READSET:END -->"

#: A short-hash token: 8 lowercase hex chars. ``20\d{6}`` is excluded so a bare date is not probed.
HASH_RE = re.compile(r"^[0-9a-f]{8}$")
DATE_LIKE_RE = re.compile(r"^20\d{6}$")
BACKTICKED_RE = re.compile(r"`([^`]+)`")
LAST_VERIFIED_RE = re.compile(r"^LAST VERIFIED:\s*(\d{4}-\d{2}-\d{2})\b", re.MULTILINE)
PLACEHOLDER_RE = re.compile(r"\b(tbd|tba|unknown|unrecorded|n/?a|todo|\?\?+)\b", re.IGNORECASE)


def _fenced(text: str, begin: str, end: str) -> str | None:
    """Return the text strictly between two fence markers, or None if absent/misordered."""
    i, j = text.find(begin), text.find(end)
    if i < 0 or j < 0 or j <= i:
        return None
    return text[i + len(begin) : j]


_ESCAPED_PIPE = "\x00PIPE\x00"


def _table_rows(block: str) -> list[list[str]]:
    """Parse pipe-table body rows out of a block, dropping header and separator lines.

    Backslash-escaped pipes (``\\|``, as in ``max\\|PF\\|``) are cell content, not delimiters.
    """
    rows: list[list[str]] = []
    for line in block.splitlines():
        s = line.strip().replace("\\|", _ESCAPED_PIPE)
        if not s.startswith("|"):
            continue
        cells = [c.strip().replace(_ESCAPED_PIPE, "|") for c in s.strip("|").split("|")]
        if all(set(c) <= set("-: ") for c in cells):  # separator
            continue
        rows.append(cells)
    return rows[1:] if rows else rows  # first surviving row is the header


def _budgeted_bytes(text: str) -> tuple[int, int]:
    """Return ``(bytes counted against MAX_BYTES, bytes exempt as section (c))``.

    Section (c) only ever grows by its own rule, so counting it against a fixed budget makes the two
    rules mutually unsatisfiable -- see MAX_BYTES. If the fence cannot be found the exemption is
    ZERO, so a malformed file is budgeted strictly rather than leniently.
    """
    block = _fenced(text, BLOCKLIST_BEGIN, BLOCKLIST_END)
    exempt = len(block.encode("utf-8")) if block is not None else 0
    return len(text.encode("utf-8")) - exempt, exempt


def _commit_exists(sha: str) -> bool:
    """True if ``sha`` resolves to a commit object in this repository."""
    return (
        subprocess.run(
            ["git", "cat-file", "-e", f"{sha}^{{commit}}"],
            cwd=REPO_ROOT,
            capture_output=True,
        ).returncode
        == 0
    )


def _is_shallow() -> bool:
    """True if this is a shallow clone, where historical commit objects are absent.

    A shallow checkout (``actions/checkout`` defaults to ``fetch-depth: 1``) cannot resolve any
    historical hash, so the commit-existence check would fail on every row for an environment
    reason rather than a content reason. It is SKIPPED and announced instead — run this gate on a
    full clone (``fetch-depth: 0``) for the check to be enforcing.
    """
    r = subprocess.run(
        ["git", "rev-parse", "--is-shallow-repository"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    return r.stdout.strip() == "true"


def _hash_tokens(text: str) -> list[str]:
    """All backticked short-hash-looking tokens in ``text`` (dates excluded)."""
    out = []
    for tok in BACKTICKED_RE.findall(text):
        t = tok.strip()
        if HASH_RE.match(t) and not DATE_LIKE_RE.match(t):
            out.append(t)
    return out


def check(path: Path) -> list[str]:
    """Run every assertion against ``path``. Returns a list of findings (empty = PASS)."""
    findings: list[str] = []
    if not path.is_file():
        return [f"MISSING_FILE {path} does not exist"]

    text = path.read_text(encoding="utf-8")
    n_lines = len(text.splitlines())
    n_bytes = len(text.encode("utf-8"))
    n_budgeted, n_exempt = _budgeted_bytes(text)

    # --- 3. shape -------------------------------------------------------------------------
    if n_budgeted > MAX_BYTES:
        findings.append(
            f"TOO_LONG {n_budgeted:,} budgeted bytes > {MAX_BYTES:,} ({n_bytes:,} total, {n_exempt:,} of it "
            f"exempt as section (c); ~{n_bytes // 4:,} tokens, {n_lines} lines): "
            f"STATE.md is a status file, not audit doc #255. Re-wrapping does not help — the budget is "
            f"what a context window pays (proxy assertion)"
        )
    for sec in REQUIRED_SECTIONS:
        if text.count(sec) != 1:
            findings.append(f"SECTION_MISSING_OR_DUPLICATED expected exactly one '{sec}' heading")

    m = LAST_VERIFIED_RE.search(text)
    if not m:
        findings.append("NO_LAST_VERIFIED expected a line 'LAST VERIFIED: YYYY-MM-DD'")

    # --- 1. tier-(a) build commits -------------------------------------------------------
    tier_a = _fenced(text, TIER_A_BEGIN, TIER_A_END)
    if tier_a is None:
        findings.append(f"NO_TIER_A_FENCE expected {TIER_A_BEGIN} ... {TIER_A_END}")
    else:
        rows = _table_rows(tier_a)
        if not rows:
            findings.append("EMPTY_TIER_A no result rows inside the tier-(a) fence")
        for k, cells in enumerate(rows, start=1):
            if len(cells) != 4:
                findings.append(
                    f"TIER_A_ROW_SHAPE row {k}: expected 4 columns "
                    f"(result | population/device | build commit | artifact), got {len(cells)}"
                )
                continue
            result, popdev, build, artifact = cells
            for name, cell in (("result", result), ("population/device", popdev), ("artifact", artifact)):
                if not cell:
                    findings.append(f"TIER_A_EMPTY_CELL row {k}: '{name}' is empty")
            if not build:
                findings.append(f"MISSING_BUILD_COMMIT row {k}: build-commit cell is empty")
                continue
            if PLACEHOLDER_RE.search(build):
                findings.append(
                    f"PLACEHOLDER_BUILD_COMMIT row {k}: build cell is a placeholder ({build!r}); "
                    f"a result whose build cannot be established belongs in section (c)"
                )
            shas = _hash_tokens(build)
            if not shas:
                findings.append(
                    f"NO_BUILD_COMMIT_HASH row {k}: build cell names no 8-hex commit "
                    f"(cell = {build!r})"
                )

    # --- git-resolvability of every hash quoted anywhere -------------------------------
    if _is_shallow():
        print(
            "[state] NOTE — shallow clone: the commit-existence check is SKIPPED (not enforcing). "
            "Use fetch-depth: 0 to enforce it."
        )
    else:
        for sha in sorted(set(_hash_tokens(text))):
            if not _commit_exists(sha):
                findings.append(f"COMMIT_NOT_FOUND `{sha}` does not resolve to a commit in this repo")

    # --- 2. read-set paths ---------------------------------------------------------------
    readset = _fenced(text, READSET_BEGIN, READSET_END)
    if readset is None:
        findings.append(f"NO_READSET_FENCE expected {READSET_BEGIN} ... {READSET_END}")
    else:
        items = [ln for ln in readset.splitlines() if re.match(r"^\s*\d+\.\s", ln)]
        if len(items) < 8:
            findings.append(
                f"READSET_TOO_SMALL {len(items)} numbered items; the audit's minimal correct entry "
                f"point is 8"
            )
        for tok in BACKTICKED_RE.findall(readset):
            t = tok.strip()
            if "/" not in t and not t.endswith((".md", ".yaml", ".yml", ".py")):
                continue
            if HASH_RE.match(t):
                continue
            if not (REPO_ROOT / t).exists():
                findings.append(f"READSET_PATH_MISSING `{t}` is in the read set but not on disk")

    findings += _check_blocklist_refs()
    return findings


#: Files that cite the (c) blocklist by number.  A citation is a promise that the reader can look the
#: entry up; a number that resolves to a DIFFERENT entry is worse than no citation, because it sends
#: them to a real paragraph about something else.
#: A citation is `(c) 14`.  A list marker is `(a) ... (b) ... (c) 0.21 pN` or `(c) 100%` — the digits there
#: are a VALUE, so they are followed by a decimal point or a percent sign.  Excluding on that is a property
#: of the text, not a numeric range: the first draft used `n > max(defined) + 20`, an invented threshold
#: that silently swallowed its own negative control at (c) 44.
_CREF_RE = re.compile(r"\(c\)\s*(\d{1,3})(?![\d.%])")


def _check_blocklist_refs() -> list[str]:
    """Every `(c) N` reference in the repo must resolve to an entry that exists in STATE_NONQUOTABLE.md.

    WHY.  Found 2026-07-29: a Notion dev-log entry cited the retired `observe` numbers as "(c) 15" when
    they are (c) 14 — and (c) 15 is a different retraction entirely (the sf_motor cost multiplier).  A
    reader following it lands on a real, confident paragraph about the wrong thing, which is a worse
    failure than a dangling link.  Notion is out of reach of any gate; the repo need not be.

    Deliberately NOT flagged: numbers that resolve, and list markers like "(a) ... (c) 0.21 pN" that the
    regex catches — those are excluded by requiring the entry to be absent from a populated entry set, so
    a malformed STATE_NONQUOTABLE cannot make every reference look wrong.
    """
    nq = REPO_ROOT / "STATE_NONQUOTABLE.md"
    if not nq.is_file():
        return []
    defined = {int(m.group(1)) for m in re.finditer(r"(?m)^\*\*(\d+)\.\*\*", nq.read_text(encoding="utf-8"))}
    if not defined:                      # nothing parsed: report that rather than every reference
        return ["BLOCKLIST_UNPARSEABLE STATE_NONQUOTABLE.md yielded no numbered entries"]
    # SCOPE, and why it is a scope rather than a cleverer regex.  `(c) 14` is a citation; `(c) 100x100 um`
    # and `(c) 0.21 pN` are list markers whose number is a VALUE.  Two attempts to separate them from the
    # text alone failed — a numeric range swallowed its own negative control, and "is the word STATE
    # nearby" rejected 24 of 33 real citations.  What actually distinguishes them is the DOCUMENT: the
    # blocklist is cited by the files that govern what may be quoted, never by paper notes under
    # `corpus/`/`lenses/`, which use (a)/(b)/(c) to enumerate a study's setups.
    governed = [REPO_ROOT / n for n in
                ("STATE.md", "ROADMAP.md", "CLAUDE.md", "AGENTS.md", "STRUCTURE.md")]
    governed += sorted((REPO_ROOT / "aleph" / "docs").glob("STATE_LOG.md"))
    governed += sorted((REPO_ROOT / "aleph" / "outputs").rglob("REPORT.md"))
    out: list[str] = []
    for md in governed:
        if not md.is_file():
            continue
        try:
            text = md.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for m in _CREF_RE.finditer(text):
            n = int(m.group(1))
            if n in defined:
                continue
            rel = md.relative_to(REPO_ROOT)
            out.append(f"BLOCKLIST_REF_UNDEFINED {rel} cites (c) {n}, which STATE_NONQUOTABLE.md does not define")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Structural gate for the root STATE.md.")
    ap.add_argument("--file", type=Path, default=DEFAULT_STATE, help="STATE.md path to check")
    args = ap.parse_args()

    findings = check(args.file)
    if findings:
        print(f"[state] FAIL — {len(findings)} finding(s) in {args.file}")
        for f in findings:
            print(f"  {f}")
        return 1
    text = args.file.read_text(encoding="utf-8")
    rows = _table_rows(_fenced(text, TIER_A_BEGIN, TIER_A_END) or "")
    n_budgeted, n_exempt = _budgeted_bytes(text)
    date = LAST_VERIFIED_RE.search(text).group(1)  # type: ignore[union-attr]
    print(
        f"[state] OK — {args.file.name}: {len(text.splitlines())} lines, "
        f"{n_budgeted:,}/{MAX_BYTES:,} budgeted bytes (+{n_exempt:,} exempt in (c), PI 2026-08-21), "
        f"6 sections, "
        f"{len(rows)} tier-(a) rows each with a resolving build commit, "
        f"{len(set(_hash_tokens(text)))} quoted commits all resolve, "
        f"read-set paths all present, last verified {date}."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
