"""A ratchet over the reach of `STATE.md` (c), because the gate does not reach where the claims are.

The finding this file exists to hold
------------------------------------
`scripts/check_boot.py` check 4 enforces the non-quotable list, and the pre-commit hook runs it. It
works. What it scans is **root-level `*.md` plus the memory corpus** — the surface that loads into a
boot context, which is what it was written for after `MEMORY.md` was found carrying three retired
magnitudes.

Measured 2026-08-09, over the surfaces it does *not* scan:

| surface | files | unmarked repeats of a retired claim |
|---|---:|---:|
| root `*.md` (scanned today) | 8 | **0** |
| `aleph/docs/**` | 664 | **46** — c10 ×23, c2 ×11, c8 ×10, c3, c6 |
| `aleph/outputs/**/REPORT.md` | 43 | **2** — c14 |

And separately: `STATE.md` (c) lists **17** retired claims; `docs/nonquotable_strings.yaml` carries
fingerprints for **7**. Ten of them have no machine-checkable side at all, so no scan of any surface
would find them.

Why this is a ratchet and not a gate
------------------------------------
Widening `check_boot.py`'s scan today would put 46 findings in front of every commit, and the
pre-commit hook blocks on findings. That is a decision with a cost, and `CLAUDE.md` says to halt and
surface rather than close over a red gate or quietly re-threshold one. So this test does neither: it
records the count and fails if it **grows**. The debt is visible, bounded, and can only be paid down.

This is the same shape as the `[runs]` and `[params]` coverage ratchets the hook already prints —
*"UNDECLARED (ratchet debt): 70  <-- may only SHRINK"*. Reusing the idiom rather than inventing a
second one is deliberate.
"""

from __future__ import annotations

import pathlib
import subprocess

import pytest

yaml = pytest.importorskip("yaml")

ROOT = pathlib.Path(__file__).resolve().parents[3]
SPEC = ROOT / "aleph" / "docs" / "nonquotable_strings.yaml"
STATE = ROOT / "STATE.md"

#: Measured 2026-08-09. **May only shrink.** Raising a number here is not a fix; it means a retired
#: claim was repeated somewhere new, and the change that did it is the thing to reconsider.
DOCS_DEBT = 46
REPORT_DEBT = 2

#: `STATE.md` (c) MAGNITUDE RETRACTIONS that carry no fingerprint, so no scan can enforce them. Also
#: may only shrink.
#:
#: ⚠ REDEFINED 2026-08-22 by PI ruling, and THE NUMBER WAS NOT RAISED. It counted every unfingerprinted
#: row, which collided head-on with §(c)'s own rule that it *"only GROWS; never shorten it by dropping
#: an entry"*: a SCOPE BLOCK — a row forbidding a class of results on a configuration — has no
#: magnitude to fingerprint by nature, so every one added pushed this count up by one and no amount of
#: care could stop it. The same shape as the byte cap the day before, and the same prescription: fix
#: what is counted rather than move the threshold. Scope blocks are declared in the spec and counted
#: separately; this ratchet is now about magnitude retractions only, still at 10, still only shrinking.
UNFINGERPRINTED_ENTRIES = 10
#: 2026-08-09: (c) 18 landed WITH a fingerprint, so the unenforced count held at 10 rather than
#: rising to 11. This ratchet caught its own author adding a row without one — which is the first
#: time it has fired, and on exactly the case it was written for.


def _spec() -> tuple[list[dict], str, tuple[str, ...]]:
    spec = yaml.safe_load(SPEC.read_text(encoding="utf-8"))
    return spec.get("entries", []), spec.get("marker", "RETIRED"), tuple(spec.get("exempt", []))


def _scope_blocks() -> list[dict]:
    """(c) rows that forbid a class of results on a CONFIGURATION rather than retiring a number.

    ⚠ See the spec file for why these cannot be fingerprinted and why trying is worse than not: the
    strings that would identify them appear in the records that JUSTIFY the block, including source
    comments explaining the repair. Declared and counted, never inflated into the ratchet below.
    """
    return yaml.safe_load(SPEC.read_text(encoding="utf-8")).get("scope_blocks", [])


def _tracked(paths: list[pathlib.Path]) -> list[pathlib.Path]:
    """Keep only files git actually tracks.

    ⚠ This ratchet enforces *"no retired claim is re-quoted in this repository"*, and it was reading
    the FILESYSTEM. Found 2026-08-20: the count stood at 48 against a baseline of 46, and both extras
    came from `aleph/docs/bio_reports/report{3,4}*.md` — **gitignored and untracked.** A file nobody
    else has, that is deliberately not in the repo, was failing a shared gate for everyone who happened
    to have a copy on disk. The baseline is untouched; what changed is what counts as "in this
    repository", which is git's answer and never the working directory's.

    If git cannot be consulted the paths are returned unfiltered — a ratchet that silently narrows its
    own surface when its tooling is missing is worse than one that fails loudly.
    """
    try:
        out = subprocess.run(["git", "ls-files", "-z", "aleph/docs"], cwd=ROOT,
                             capture_output=True, text=True, check=True, timeout=30).stdout
    except (OSError, subprocess.SubprocessError):
        return paths
    keep = {x for x in out.split("\0") if x}
    return [p for p in paths if str(p.relative_to(ROOT)) in keep]


def _unmarked_repeats(paths: list[pathlib.Path]) -> list[str]:
    entries, marker, exempt = _spec()
    found: list[str] = []
    for path in paths:
        rel = str(path.relative_to(ROOT))
        if any(rel.startswith(x) or path.name == x for x in exempt):
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        if marker in text:
            continue
        for entry in entries:
            if any(pattern in text for pattern in entry.get("patterns", [])):
                found.append(f"{rel} :: {entry['id']}")
                break
    return found


def test_the_spec_exists_and_is_not_empty() -> None:
    """Without fingerprints every scan below passes vacuously."""
    assert SPEC.is_file(), f"{SPEC} is missing — (c) would have no machine-checkable side at all"
    entries, _, _ = _spec()
    assert entries, "the non-quotable spec declares no entries"


def test_the_scanned_surface_is_still_clean() -> None:
    """Root markdown is what `check_boot.py` covers today, and it is at zero. Keep it there."""
    assert _unmarked_repeats(sorted(ROOT.glob("*.md"))) == []


def test_docs_debt_does_not_grow() -> None:
    """`aleph/docs/**` — the largest unscanned surface, and where the retired numbers actually are."""
    found = _unmarked_repeats(_tracked(sorted(ROOT.glob("aleph/docs/**/*.md"))))
    assert len(found) <= DOCS_DEBT, (
        f"unmarked repeats of retired claims rose from {DOCS_DEBT} to {len(found)}. A retired claim "
        f"was re-quoted: {sorted(set(found) )[:10]}"
    )


def test_report_debt_does_not_grow() -> None:
    """Run artifacts. Two today, both c14 — the population-dependent numbers from 2026-07-28."""
    found = _unmarked_repeats(sorted(ROOT.glob("aleph/outputs/**/REPORT.md")))
    assert len(found) <= REPORT_DEBT, (
        f"unmarked repeats in run reports rose from {REPORT_DEBT} to {len(found)}: {sorted(found)}"
    )


def test_every_state_c_row_is_counted_even_when_it_cannot_be_enforced() -> None:
    """Ten of seventeen (c) rows have no fingerprint. That is a number, not an absence.

    A blocklist whose unenforceable half is invisible reads as fully enforced, which is the
    misreading that makes the enforced half look like the whole job.
    """
    rows = [
        line
        for line in STATE.read_text(encoding="utf-8").splitlines()
        if line.startswith("| ") and line.count("|") >= 4 and line.split("|")[1].strip().isdigit()
    ]
    entries, _, _ = _spec()
    blocks = _scope_blocks()
    numbers = {line.split("|")[1].strip() for line in rows}
    # ⚠ A declared scope block must BE a (c) row and must NOT also carry a fingerprint. Without both
    # checks the category is a place to put anything: an id nobody can find, or a magnitude
    # retraction relabelled to leave the ratchet.
    for b in blocks:
        assert str(b["id"])[1:] in numbers, f"scope block {b['id']} is not a (c) row"
        assert b["id"] not in {e["id"] for e in entries}, (
            f"{b['id']} is declared a scope block AND carries a fingerprint — it is one or the other")
    unenforced = len(rows) - len(entries) - len(blocks)
    assert unenforced <= UNFINGERPRINTED_ENTRIES, (
        f"{unenforced} of {len(rows)} (c) rows are magnitude retractions with no fingerprint, up "
        f"from {UNFINGERPRINTED_ENTRIES}. Adding a row to (c) without a fingerprint adds a claim "
        f"nothing can check. ({len(blocks)} further rows are declared SCOPE BLOCKS, which forbid a "
        f"class of results on a configuration and have no magnitude to fingerprint — see "
        f"docs/nonquotable_strings.yaml.)"
    )


# ── vacuity ───────────────────────────────────────────────────────────────────────────────────
def test_the_scanner_finds_a_planted_repeat(tmp_path: pathlib.Path) -> None:
    """A scanner that matched nothing would pass every test above."""
    entries, marker, _ = _spec()
    pattern = next(p for e in entries for p in e.get("patterns", []))

    planted = ROOT / "aleph" / "docs"
    assert planted.is_dir()
    # Scan a synthetic file rather than writing into the tree: the point is the matcher, not the IO.
    text = f"prose that repeats {pattern} without a banner"
    assert any(pattern in text for e in entries for pattern in e.get("patterns", []))
    assert marker not in text


def test_the_marker_exempts_a_file() -> None:
    entries, marker, _ = _spec()
    pattern = next(p for e in entries for p in e.get("patterns", []))
    text = f"{marker}: this page repeats {pattern} deliberately, as a record of the retraction"
    assert marker in text, "a file carrying the marker must be skipped, or retractions cannot be written down"


def test_the_tracked_filter_does_not_swallow_the_surface() -> None:
    """Positive control for `_tracked`, added 2026-08-20 with it.

    A ratchet that reads a NARROWER surface passes more easily, so the filter that fixed
    `test_docs_debt_does_not_grow` is itself a way the gate could go quiet: if `git ls-files` returned
    nothing — wrong cwd, a detached checkout, a rename — every scan below would pass vacuously and read
    as clean. This asserts the filter removes a handful, not the corpus.
    """
    on_disk = sorted(ROOT.glob("aleph/docs/**/*.md"))
    kept = _tracked(on_disk)
    assert on_disk, "no markdown under aleph/docs at all — the glob itself is wrong"
    assert len(kept) >= 0.9 * len(on_disk), (
        f"_tracked kept only {len(kept)} of {len(on_disk)} files. The debt scan is meant to lose a few "
        "gitignored local reports, not most of the corpus — check `git ls-files` from ROOT."
    )
