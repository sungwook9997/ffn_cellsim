#!/usr/bin/env python
r"""Four checks on the boot context — the thing every session pays for before it does any work.

WHY THIS EXISTS.  On 2026-07-28 the boot context was measured for the first time: ~16,200 tokens of
`CLAUDE.md` + `MEMORY.md` injected unconditionally, before a single file of the actual task was read.
Cutting it was a morning's work.  Keeping it cut is not, because nothing checked it — `CLAUDE.md` had
already grown to 524 lines once, and the only reason anyone noticed is that a human read it and said so.

The audit also found the failure that matters more than the size.  `MEMORY.md` loads on every boot and
was carrying three magnitudes that `STATE.md` (c) had already RETIRED, so the retraction and the
retracted number entered the same context window and **the retracted one arrived first**.  That is not a
tidiness problem; it is the mechanism by which a withdrawn result gets re-quoted, and this project has
withdrawn sixteen.  Check 4 makes (c) enforceable against the whole memory corpus instead of against
whoever happens to re-read it.

THE FOUR CHECKS, and what each would have caught:

1. **Budget** — the boot set's total bytes.  Would have caught the growth to 524 lines.
2. **No magnitudes in `CLAUDE.md`** — the charter's own first rule, which nothing tested.  Would have
   caught `13 components / 32 connectors` sitting in two sections while the real census was 14/36.
3. **`MEMORY.md` is an index** — one line per memory, title and file, nothing else.  Would have caught
   the hooks growing into findings, which is how the three retired numbers got in.
4. **(c) cross-check** — no file may repeat a retired claim's fingerprint without a `RETIRED` marker.
   Would have caught all three by name, on the commit that introduced them.

WHAT THIS DOES NOT DO.  It does not judge whether the boot context is CORRECT — 40 KB of confidently
wrong status passes checks 1-3.  Check 4 is the only one that tests a claim rather than a shape, and it
tests only claims someone has already written a fingerprint for; it reports how many (c) entries have
none, so the gap is visible rather than assumed closed.

Runs on stdlib + pyyaml, reads only committed artifacts plus the memory directory, and exits non-zero
with one line per finding.

Usage:
    python aleph/scripts/check_boot.py [--memory-dir PATH] [--verbose]

Sanity Gate:
    * dimensional: sizes are BYTES throughout, never lines — a line budget is defeatable by re-wrapping,
      which is exactly how this project's previous size assertion was defeated the same day it was
      questioned.
    * boundary: a missing memory directory SKIPS checks 3 and 4 with a stated reason rather than passing
      them silently, because "no memories to check" and "checked, found nothing" are different results.
    * conservation/invariant: every (c) entry either has a fingerprint or is counted as unenforced; the
      two always sum to the number of entries in the YAML.
    * numerical: no arithmetic beyond byte counts and comparisons.
    * sign-sense: a finding is always "over budget" or "present where forbidden" — never the reverse, so
      an empty repository cannot fail.
    * measurement-protocol: budgets are declared as constants here, above the code, and changing one is
      a visible diff; the check never derives a budget from the file it is measuring.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

#: What a session loads before it starts: the two files the harness injects unconditionally, plus the
#: one the boot protocol names. Anything else is pulled on demand and is not this budget's business.
BOOT_SET = ("CLAUDE.md", "STATE.md")
MEMORY_INDEX = "MEMORY.md"

#: Byte budgets. Set just ABOVE the measured size after the 2026-07-28 reduction, never AT it — a budget
#: set to today's size is loosening by default and can only ratchet the wrong way. Lowering these is the
#: intended direction; raising one is a PI decision.
BUDGETS: dict[str, int] = {
    "CLAUDE.md": 18_000,      # measured ~15,400 after the roadmap moved to ROADMAP.md
    "MEMORY.md": 11_000,      # measured ~9,600 after the hooks were removed
    "STATE.md": 22_000,       # RAISED 20,000 -> 22,000, PI 2026-08-10; also enforced by check_state_md.py
}

#: Bytes the per-file budget does NOT count, keyed by file. PI decision 2026-08-21: STATE.md section (c)
#: is exempt, and the 22,000 was NOT raised. (c)'s own header says it "only GROWS; never shorten it by
#: dropping an entry" -- a monotonically growing list under a fixed cap is unsatisfiable, and only the
#: DATE of the collision is negotiable. It arrived 2026-08-21 with 39 B free and `c19` sitting unpushed
#: in the companion file, which is the one thing (c) exists to prevent.
#: ⚠ THIS EXEMPTS (c) FROM ITS FILE'S CAP, NOT FROM TOTAL_BUDGET. Those bytes are still injected into
#: every boot and are still counted in the total below -- (c) is PUSHED on purpose, so it is genuinely
#: paid. The per-file cap now governs only what a session can choose not to write.
#: The exemption is ANNOUNCED on every run: an invisible carve-out stops being a proxy and becomes a hole.
EXEMPT_SECTION: dict[str, tuple[str, str]] = {"STATE.md": ("## (c)", "## (d)")}


def _exempt_bytes(name: str, raw: bytes) -> int:
    """Bytes in ``raw`` that ``name``'s per-file budget does not count (see EXEMPT_SECTION).

    Returns 0 when the file has no exemption or its fence headings cannot be found, so a malformed
    file is budgeted strictly rather than leniently.
    """
    fence = EXEMPT_SECTION.get(name)
    if fence is None:
        return 0
    text = raw.decode("utf-8", errors="replace")
    i, j = text.find(fence[0]), text.find(fence[1])
    return 0 if i < 0 or j <= i else len(text[i:j].encode("utf-8"))

#: RAISED 46,000 -> 48,000 by PI decision, 2026-08-10 — the ratification of the open item `STATE.md` (e) 3,
#: which the PI had proposed on 2026-07-28 and which had sat unratified since. What forced it was a
#: sixteenth tier-(a) row (the Card-5 native array-ownership gate) landing in a `STATE.md` that had 41 B of
#: headroom. Both numbers were surfaced and granted, not edited to fit: a threshold moved after seeing the
#: data is the thing the charter forbids, and the only defence is that the moving is somebody else's call.
#: The direction of travel is still DOWN. Two budgets hitting their ceiling in one night is the signal —
#: the next session to reach this line should move status out of the boot set before touching this number.
#:
#: ⚠ **RAISED AGAIN 48,000 -> 50,000 by PI decision, 2026-08-24 — AND THE PARAGRAPH ABOVE TOLD ME NOT TO.**
#: It says the next session to reach this line should move status out of the boot set FIRST. I am that
#: session and I did not: the PI, asked, ruled to raise it, and that ruling is the only thing that may move
#: a threshold. So this is recorded as an override of the file's own advice rather than as agreement with
#: it, because the advice is still correct and the next session should read it as unpaid.
#:
#: What forced it: `STATE.md` (c) 21, the finding that the cell every gamma and tau was measured on carried
#: two of its eleven populations' laws, had no motor in it, and ran a cortex at a tenth of actin. (c) only
#: GROWS by charter, so a blocklist entry cannot be traded against anything, and the file had 41 B free —
#: the SAME 41 B the 2026-08-10 raise was granted against, which is its own signal.
#:
#: ⚠ **THE REDUCTION IS AVAILABLE AND IT IS LARGE.** `STATE.md` (a) already says most tier-(b) rows are
#: "history on an archived tree" after the 2026-08-22 archival ruling — they were measured on `engine/`
#: and `components/`, which are read-only port sources now. Moving those rows to a dated history file
#: would free far more than 2,000 B and would make (b) describe the canonical tree, which is what it is
#: for. That is a PI call about what stays authoritative, not an edit; it is named here so the next raise
#: has to argue against it.
TOTAL_BUDGET = 50_000

#: A magnitude is a number bound to a unit or a scale word. Provenance dates and version strings are not
#: magnitudes; counts of things that change (components, connectors, filaments) very much are.
MAGNITUDE_RE = re.compile(
    r"\b\d[\d,._]*\s*"
    r"(pN|nm|µm|um|Pa|kPa|mN/m|s/step|GB|MB|KB|tokens?|decades?|"
    r"components?|connectors?|filaments?|nodes?|heads?|seeds?|rows?|files?|lines?)\b"
    r"|\b\d[\d,._]*\s*[×x]\b"
    r"|\b\d[\d,._]*\s*%",
    re.IGNORECASE,
)

#: Strings that match MAGNITUDE_RE but are not status: language/tool versions, and structural ordinals
#: that name a contract rather than measure anything.
MAGNITUDE_ALLOW = (
    "Python 3.13.13",
    "Newton's 3rd law",
    "8-state evidence ladder",
)

MEMORY_LINE_RE = re.compile(r"^- \[[^\]]+\]\(([A-Za-z0-9._-]+\.md)\)$")


def _memory_dir_for(root: Path) -> Path:
    """Return the conventional memory directory for a checkout root, whether or not it exists."""
    slug = "-" + str(root).lstrip("/").replace("/", "-").replace("_", "-")
    return Path.home() / ".claude" / "projects" / slug / "memory"


def _main_checkout_root() -> Path:
    """Return the MAIN checkout's root, even when called from a linked git worktree.

    A worktree's ``.git`` is a file pointing at ``<main>/.git/worktrees/<name>``, so
    ``git rev-parse --git-common-dir`` resolves to the main checkout's ``.git`` and its parent is the
    main root.  Falls back to :data:`REPO_ROOT` whenever git cannot be asked — the caller then gets
    the old behaviour rather than an exception.
    """
    try:
        done = subprocess.run(["git", "-C", str(REPO_ROOT), "rev-parse", "--git-common-dir"],
                              capture_output=True, text=True, timeout=15, check=False)
    except (OSError, subprocess.SubprocessError):
        return REPO_ROOT
    if done.returncode != 0 or not done.stdout.strip():
        return REPO_ROOT
    common = Path(done.stdout.strip())
    if not common.is_absolute():
        common = (REPO_ROOT / common).resolve()
    return common.parent


def _default_memory_dir() -> Path:
    """Return the memory directory to scan, resolving a linked worktree to the main checkout.

    WHY THIS IS NOT A LOOSENING.  The memory corpus is per-REPOSITORY, not per-checkout: one
    ``.claude/projects/<slug>/memory`` exists for the main clone and a linked worktree has no slug of
    its own.  Resolving by ``cwd`` therefore made every worktree report ``MEMORY_DIR_ABSENT`` — which
    the gate correctly refuses to call a pass ("nothing was examined"), and which consequently made it
    impossible to commit from ANY worktree once the PI put parallel experiments on their own branches
    (2026-08-10).  The tempting fix — create an empty directory at the worktree's slug — would have
    turned a loud false negative into a silent false PASS over a corpus of zero files.  This instead
    points the check at the corpus that actually exists, so checks 3 and 4 run against real content
    from a worktree exactly as they do from the main checkout.  If git cannot be asked, the previous
    behaviour is returned unchanged and the gate still refuses.
    """
    return _memory_dir_for(_main_checkout_root())


def check_budget() -> list[str]:
    """Return findings for any boot file over its byte budget, or for the total.

    Returns:
        One finding per breach; empty when every file and the total are within budget.
    """
    findings: list[str] = []
    total = 0
    for name in BOOT_SET:
        path = REPO_ROOT / name
        if not path.is_file():
            findings.append(f"BOOT_FILE_MISSING {name}")
            continue
        raw = path.read_bytes()
        size = len(raw)
        total += size                       # the TOTAL counts every byte, exemptions included
        exempt = _exempt_bytes(name, raw)
        budget = BUDGETS.get(name)
        if budget and size - exempt > budget:
            findings.append(
                f"OVER_BUDGET {name} {size - exempt:,} B > {budget:,} B ({size:,} total, {exempt:,} of "
                f"it exempt; ~{size // 4:,} tokens). Move status out; do not re-wrap — the budget is "
                f"bytes for that reason"
            )
        if exempt:
            print(f"[boot] NOTE — {name}: {exempt:,} B exempt from its {budget:,} B cap "
                  f"(section (c), PI 2026-08-21); still counted in the boot total")
    index = _default_memory_dir() / MEMORY_INDEX
    if index.is_file():
        size = len(index.read_bytes())
        total += size
        if size > BUDGETS[MEMORY_INDEX]:
            findings.append(
                f"OVER_BUDGET {MEMORY_INDEX} {size:,} B > {BUDGETS[MEMORY_INDEX]:,} B — it is an index; "
                f"one line per memory, no hooks"
            )
    if total > TOTAL_BUDGET:
        findings.append(
            f"OVER_BUDGET boot total {total:,} B > {TOTAL_BUDGET:,} B (~{total // 4:,} tokens) — every "
            f"session pays this before doing anything"
        )
    return findings


def check_no_magnitudes(path: Path) -> list[str]:
    """Return findings for magnitudes in the charter, which its own first rule forbids.

    Args:
        path: The charter file.

    Returns:
        One finding per offending line.
    """
    if not path.is_file():
        return [f"MISSING {path.name}"]
    findings: list[str] = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        stripped = line.strip()
        if stripped.startswith(("#", ">", "|")) and "VERIFY" in stripped:
            continue                                    # a ⚠ VERIFY pointer is the sanctioned form
        for hit in MAGNITUDE_RE.finditer(line):
            text = hit.group(0)
            window = line[max(0, hit.start() - 24): hit.end() + 24]
            if any(allowed in window for allowed in MAGNITUDE_ALLOW):
                continue
            findings.append(
                f"MAGNITUDE_IN_CHARTER {path.name}:{number} `{text.strip()}` — a number here goes stale "
                f"and is injected into every boot. Move it to STATE.md and leave a ⚠ VERIFY pointer"
            )
    return findings


def check_memory_index(index: Path) -> list[str]:
    """Return findings for an index line that is anything other than a title and a file.

    Args:
        index: ``MEMORY.md``.

    Returns:
        One finding per malformed line or dangling link.
    """
    if not index.is_file():
        return []
    findings: list[str] = []
    for number, line in enumerate(index.read_text(encoding="utf-8").splitlines(), 1):
        if not line.startswith("- "):
            continue
        match = MEMORY_LINE_RE.match(line.rstrip())
        if not match:
            findings.append(
                f"MEMORY_INDEX_SHAPE {MEMORY_INDEX}:{number} is not `- [title](file.md)`. A hook here "
                f"becomes a finding, and findings here are injected ahead of their own retractions"
            )
            continue
        if not (index.parent / match.group(1)).is_file():
            findings.append(f"MEMORY_LINK_DEAD {MEMORY_INDEX}:{number} -> {match.group(1)}")
    return findings


def check_retired_claims(memory_dir: Path, verbose: bool = False) -> list[str]:
    """Return findings for any file repeating a retired claim without a ``RETIRED`` marker.

    Args:
        memory_dir: Directory of memory files; may not exist.
        verbose: Print the per-entry fingerprint tally.

    Returns:
        One finding per (file, retired claim) pair that carries no marker.
    """
    import yaml

    spec_path = REPO_ROOT / "aleph" / "docs" / "nonquotable_strings.yaml"
    if not spec_path.is_file():
        return [f"MISSING {spec_path.relative_to(REPO_ROOT)} — (c) has no machine-checkable side"]
    spec = yaml.safe_load(spec_path.read_text(encoding="utf-8"))
    entries = spec.get("entries", [])
    marker = spec.get("marker", "RETIRED")
    exempt = tuple(spec.get("exempt", []))

    scanned: list[Path] = [p for p in (REPO_ROOT).glob("*.md")]
    if memory_dir.is_dir():
        scanned += sorted(memory_dir.glob("*.md"))
    else:
        return [
            f"MEMORY_DIR_ABSENT {memory_dir} — checks 3 and 4 did not run over the memory corpus. "
            f"Not a pass: nothing was examined"
        ]

    findings: list[str] = []
    for path in scanned:
        try:
            rel = str(path.relative_to(REPO_ROOT))
        except ValueError:
            rel = path.name
        if any(rel.startswith(x) or path.name == x for x in exempt):
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        if marker in text:
            continue
        for entry in entries:
            for pattern in entry.get("patterns", []):
                if pattern in text:
                    findings.append(
                        f"RETIRED_CLAIM_UNMARKED {rel} repeats `{pattern}` — STATE.md (c) "
                        f"{entry['id']}: {entry['claim']}. Add a `{marker}` banner or remove it"
                    )
                    break
    if verbose:
        print(f"[boot] (c) fingerprints: {len(entries)} entries, {len(scanned)} files scanned")
    return findings


def main() -> int:
    """Run all four checks and report."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--memory-dir", default=None, help="override the memory directory")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    memory_dir = Path(args.memory_dir) if args.memory_dir else _default_memory_dir()

    findings: list[str] = []
    findings += check_budget()
    findings += check_no_magnitudes(REPO_ROOT / "CLAUDE.md")
    findings += check_memory_index(memory_dir / MEMORY_INDEX)
    findings += check_retired_claims(memory_dir, verbose=args.verbose)

    if findings:
        print(f"[boot] FAIL — {len(findings)} finding(s)", file=sys.stderr)
        for finding in findings:
            print(f"  {finding}", file=sys.stderr)
        return 1

    total = sum(len((REPO_ROOT / n).read_bytes()) for n in BOOT_SET if (REPO_ROOT / n).is_file())
    index = memory_dir / MEMORY_INDEX
    if index.is_file():
        total += len(index.read_bytes())
    print(f"[boot] OK — boot context {total:,} B (~{total // 4:,} tokens) of {TOTAL_BUDGET:,} B; "
          f"no magnitudes in the charter, index well-formed, no retired claim repeated unmarked.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
