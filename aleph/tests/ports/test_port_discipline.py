"""Mechanical enforcement of the port ledger discipline.

The rule is about intent — write the entry before the code, stand on a control this tree owns, ship
a negative control that can actually fail. Intent does not survive a long autonomous build. This
module turns the mechanisable part into tests that fail loudly and name the offending file.

The sharp one is :func:`test_named_controls_resolve_to_real_tests`. Naming a test in a ledger entry
costs nothing, so an entry that names controls proves nothing on its own;
:func:`test_accepted_entries_name_both_controls` catches an entry that never claimed a control, and
this one catches an entry that claimed one nobody wrote. Everything else here is bookkeeping by
comparison.

Carried into this tree with the merge of 2026-08-09, deliberately as the **test** and not only the
template. `Project_Aleph` had both; a template with no test is a convention, and a convention is
what this tree already had.
"""

from __future__ import annotations

import ast
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[3]
LEDGER = ROOT / "ports" / "ledger"
#: `Project_Aleph`'s 135 records of ports **out of** this tree. Deliberately not scanned by the
#: checks below: they name controls in another repository and cite commits in another template, so
#: every failure they produce is a category error rather than a defect. `ports/inherited/README.md`
#: is the argument. What IS checked about them is that nobody edits one — see the two tests at the
#: end of this module.
INHERITED = ROOT / "ports" / "inherited"
TEMPLATE = ROOT / "ports" / "TEMPLATE.md"
TESTS = ROOT / "aleph" / "tests"

#: The fourteen headings every entry must carry. Matched by number, because an entry is allowed to
#: extend the title — `ALEPH-PORT-4001` §3 ends with "-- it does not, and that is the finding",
#: which is exactly the kind of honesty the heading is asking for.
MANDATORY_SECTIONS = tuple(range(1, 15))

RECOGNISED_STATUSES = frozenset({"PROPOSED", "AUDITED", "ACCEPTED", "REJECTED"})

PLACEHOLDERS = ("TODO", "TBD", "FIXME", "???", "<slug>", "NNN")

_SECTION = re.compile(r"^##\s+(\d{1,2})\.\s", re.MULTILINE)
#: The **whole** rest of the line, not a non-greedy cell: the template separates its four
#: alternatives with escaped pipes, so a cell-wise match would stop at the first one and read the
#: unfilled template as a valid `PROPOSED`.
_STATUS_ROW = re.compile(r"^\|\s*Status\s*\|(.*)$", re.MULTILINE)
_TEST_NAME = re.compile(r"\btest_[A-Za-z0-9_]+\b")

#: A ledger's control table names a **file** and a **function**, and `tests/units/test_band.py`
#: contains the token `test_band`. Scanning the raw text would demand a function of that name and
#: fail on a correctly written entry, so paths are removed before the function scan.
_PY_PATH = re.compile(r"\S*\.py\b")


def _named_tests(text: str) -> set[str]:
    """Test function names in a section, with `*.py` paths removed first."""
    return set(_TEST_NAME.findall(_PY_PATH.sub(" ", text)))


def _defined_test_functions() -> set[str]:
    """Every `def test_*` actually defined under `tests/`, parsed rather than grepped.

    Grepping the token was the first spelling and it produced a **false pass**: this module's own
    self-test contains the string ``tests/units/test_band.py``, which put ``test_band`` into the
    "defined" set and let a ledger entry name a function nobody had written. A mention is not a
    definition, and only the parser knows the difference.
    """
    defined: set[str] = set()
    for path in TESTS.rglob("test_*.py"):
        try:
            tree = ast.parse(path.read_text(errors="ignore"))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith(
                "test_"
            ):
                defined.add(node.name)
    return defined


def ledger_entries() -> list[pathlib.Path]:
    """Every ledger entry. The template is not one and lives outside the directory."""
    return sorted(p for p in LEDGER.glob("*.md") if p.name != "INDEX.md")


def _section_text(body: str, number: int) -> str:
    """The body of `## <number>.` up to the next `##` heading, or '' if absent.

    The heading's own title is excluded: a section whose only content is the words in its heading
    is empty, and `test_ledger_entries_are_complete` has to be able to say so.
    """
    starts = {int(m.group(1)): m for m in _SECTION.finditer(body)}
    if number not in starts:
        return ""
    line_end = body.find("\n", starts[number].end())
    start = len(body) if line_end == -1 else line_end + 1
    later = [m.start() for m in _SECTION.finditer(body) if m.start() >= start]
    return body[start : (later[0] if later else len(body))]


def _declared_status(body: str) -> str:
    """The single status an entry declares, or '' if it declares none or more than one.

    Returning '' for an ambiguous cell is the point. The template's row lists all four alternatives
    separated by escaped pipes, and an entry that was filled in by deleting nothing would otherwise
    read as whichever one happens to come first.
    """
    match = _STATUS_ROW.search(body)
    if not match:
        return ""
    named = [s for s in RECOGNISED_STATUSES if s in match.group(1)]
    return named[0] if len(named) == 1 else ""


@pytest.fixture(scope="module")
def entries() -> list[tuple[pathlib.Path, str]]:
    found = ledger_entries()
    assert found, f"no ledger entries under {LEDGER} — the discipline test would pass vacuously"
    return [(p, p.read_text()) for p in found]


def test_the_template_exists_and_declares_every_mandatory_section() -> None:
    """If the template drifts from this list, entries copied from it start out non-compliant."""
    assert TEMPLATE.is_file(), f"{TEMPLATE} is missing"
    numbers = {int(m.group(1)) for m in _SECTION.finditer(TEMPLATE.read_text())}
    assert set(MANDATORY_SECTIONS) <= numbers, sorted(set(MANDATORY_SECTIONS) - numbers)


def test_ledger_entries_are_complete(entries: list[tuple[pathlib.Path, str]]) -> None:
    """Every mandatory heading present, and none of them left empty."""
    problems: list[str] = []
    for path, body in entries:
        numbers = {int(m.group(1)) for m in _SECTION.finditer(body)}
        for missing in sorted(set(MANDATORY_SECTIONS) - numbers):
            problems.append(f"{path.name}: no section {missing}")
        for number in sorted(numbers & set(MANDATORY_SECTIONS)):
            if not _section_text(body, number).strip():
                problems.append(f"{path.name}: section {number} is empty")
    assert not problems, problems


def test_no_entry_still_carries_a_placeholder(entries: list[tuple[pathlib.Path, str]]) -> None:
    """An entry written as a formality reads exactly like the template it was copied from."""
    problems = [
        f"{path.name}: {token!r}"
        for path, body in entries
        for token in PLACEHOLDERS
        if token in body
    ]
    assert not problems, problems


def test_every_entry_declares_a_recognised_status(entries: list[tuple[pathlib.Path, str]]) -> None:
    problems = [
        f"{path.name}: status {_declared_status(body)!r}"
        for path, body in entries
        if _declared_status(body) not in RECOGNISED_STATUSES
    ]
    assert not problems, problems


def test_accepted_entries_name_both_controls(entries: list[tuple[pathlib.Path, str]]) -> None:
    """An ACCEPTED entry with no positive and no deliberately-failing negative control is a claim."""
    problems: list[str] = []
    for path, body in entries:
        if _declared_status(body) != "ACCEPTED":
            continue
        for number, kind in ((8, "positive"), (9, "negative")):
            if not _named_tests(_section_text(body, number)):
                problems.append(f"{path.name}: ACCEPTED but section {number} names no {kind} control")
    assert not problems, problems


def test_named_controls_resolve_to_real_tests(entries: list[tuple[pathlib.Path, str]]) -> None:
    """The sharp end. Naming a test costs nothing; this checks somebody wrote it.

    Applies at **every** status, not only ACCEPTED — an entry that names a control it never wrote
    is misleading while it is still PROPOSED, and by the time it reaches ACCEPTED the claim has
    been read several times.
    """
    defined = _defined_test_functions()

    problems: list[str] = []
    for path, body in entries:
        named: set[str] = set()
        for number in (8, 9):
            named.update(_named_tests(_section_text(body, number)))
        for missing in sorted(named - defined):
            problems.append(f"{path.name} names {missing}, which is defined nowhere under tests/")
    assert not problems, problems


def test_entries_record_a_source_commit(entries: list[tuple[pathlib.Path, str]]) -> None:
    """Section 2 must carry a 40-character sha. A port citing a branch cites a moving target."""
    problems = [
        f"{path.name}: section 2 has no 40-char sha"
        for path, body in entries
        if not re.search(r"\b[0-9a-f]{40}\b", _section_text(body, 2))
    ]
    assert not problems, problems


# ── self-tests: a scanner that silently stops matching would pass everything ───────────────────
def test_the_section_scanner_actually_matches() -> None:
    sample = "## 1. API\nbody one\n\n## 2. Source identity\nbody two\n"
    assert {int(m.group(1)) for m in _SECTION.finditer(sample)} == {1, 2}
    assert _section_text(sample, 1).strip() == "body one"
    assert _section_text(sample, 2).strip() == "body two"
    assert _section_text(sample, 3) == ""


def test_the_status_scanner_actually_matches() -> None:
    assert _declared_status("| Status | `ACCEPTED` |") == "ACCEPTED"
    assert _declared_status("| Status | `PROPOSED` |") == "PROPOSED"
    assert _declared_status("no status row here") == ""
    # The template's alternation must NOT read as a recognised status, or the template would pass.
    assert _declared_status("| Status | `PROPOSED` \\| `AUDITED` |") not in RECOGNISED_STATUSES


def test_a_mention_is_not_a_definition() -> None:
    """Regression: the token scan this replaced counted a filename in a docstring as a definition.

    `test_band` appears in this module as part of a path in another self-test. If the resolver ever
    goes back to grepping, that mention re-enters the "defined" set and a ledger entry can name a
    function nobody wrote.
    """
    defined = _defined_test_functions()
    assert "test_a_mention_is_not_a_definition" in defined
    assert "test_band" not in defined, "a path fragment is being counted as a test definition"


def test_the_inherited_records_are_present_and_separate() -> None:
    """They travel with the mechanism (B5) but not into the active ledger.

    A record of a port into another repository is evidence about why that tree's contracts exist. It
    is not a port into this one, and scoring it against this tree's discipline says something false
    about it.
    """
    inherited = sorted(INHERITED.glob("ALEPH-PORT-*.md"))
    assert len(inherited) == 135, (
        f"{len(inherited)} inherited records; B5 carried 135. The number was stated as 136 until a "
        "count of `ls *.md` was noticed to include INDEX.md"
    )
    assert (INHERITED / "README.md").is_file(), "the inherited set needs its own rules stated"

    # And none of them has leaked into the directory the discipline test governs.
    active = {p.name for p in LEDGER.glob("ALEPH-PORT-*.md")}
    leaked = active & {p.name for p in inherited}
    assert not leaked, f"inherited records in the active ledger: {sorted(leaked)}"


def test_the_new_direction_uses_its_own_id_band() -> None:
    """Inherited IDs stop below 4000; ports into this tree start at 4001, so the ID says which way."""
    import re as _re

    def _ids(directory: pathlib.Path) -> list[int]:
        out = []
        for path in directory.glob("ALEPH-PORT-*.md"):
            match = _re.match(r"ALEPH-PORT-(\d+)", path.name)
            if match:
                out.append(int(match.group(1)))
        return out

    assert max(_ids(INHERITED)) < 4000, "an inherited record is using the new band"
    assert min(_ids(LEDGER)) >= 4001, "a port into this tree is using an inherited ID"


def test_the_control_scanner_actually_matches() -> None:
    assert _named_tests("`…::test_guard_returns_the_value_untouched` asserts") == {
        "test_guard_returns_the_value_untouched"
    }
    assert _named_tests("no controls named here") == set()
    # The false positive this scanner was written against: a path is not a function.
    assert _named_tests("`aleph/tests/units/test_band.py::test_membership_matches`") == {
        "test_membership_matches"
    }
