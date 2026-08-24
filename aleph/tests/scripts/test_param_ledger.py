"""The parameter sweep must not narrow itself — that is the only way it can lie.

⚠ **WHY THIS EXISTS.** On 2026-08-24 the parameter inventory was assembled by hand and came to 79 new
axes. An AST sweep of `laws/` alone found 185, missing whole modules: `fa_maturation` (the talin chain),
`ecm_library`, `hand_kmc`, `polarization_activegel`, `polymerization_warp` and `piezo`. A hand count of
declarations is not slightly wrong, it is wrong by more than it finds, and it is wrong again next time.

⚠ **AND THE SWEEP ITSELF HAD THE SAME DEFECT ON ITS FIRST RUN.** `_scan_module` handled `ast.Assign`
and not `ast.AnnAssign`, so every module-level ANNOTATED constant (`X: float = 1.0`) was invisible —
neither counted nor listed as excluded. It appeared in no total and in no exclusion, which is the exact
shape the exclusion list exists to prevent. Its own `_demo` caught it on `turgor_pi0.LEDGER_PATH`.
`test_the_annotated_module_constant_form_is_covered` is the ratchet on that specific hole.
"""

from __future__ import annotations

import ast

from aleph.scripts.param_ledger import ROOTS, Scan, _scan_module, sweep


def test_rows_and_exclusions_partition_the_scan() -> None:
    """Nothing the sweep names may fall between the two lists.

    A declaration that is in neither is invisible in exactly the way a hand count is, and it inflates
    confidence rather than the total — the reader sees a complete-looking ledger.
    """
    scan = sweep()
    assert scan.files > 50, scan.files
    assert scan.rows, "swept the tree and found no declaration — that is a finding, not a pass"
    named = {(r.module, r.name) for r in scan.rows}
    excluded = {(m, n) for m, n, _ in scan.excluded}
    assert not (named & excluded), sorted(named & excluded)[:5]
    assert all(reason.strip() for _m, _n, reason in scan.excluded), "an exclusion with no reason"


def test_the_annotated_module_constant_form_is_covered(tmp_path) -> None:
    """`X: float = 1.0` at module level must reach the ledger. The ratchet on the sweep's own bug.

    ⚠ Written against a SYNTHETIC module rather than a real one, so it keeps testing the parser after
    whichever real constant it would otherwise have pinned is renamed or deleted.
    """
    source = tmp_path / "probe.py"
    source.write_text(
        "PLAIN = 2.5\n"
        "ANNOTATED: float = 3.5\n"
        "ANNOTATED_INT: int = 7\n"
        "LEDGER_PATH: str = 'x/y'\n"
        "NOT_A_NUMBER = ('a', 'b')\n",
        encoding="utf-8",
    )
    scan = Scan()
    _scan_module(source, "probe.py", scan)
    found = {r.name: r.value for r in scan.rows}
    assert found == {"PLAIN": 2.5, "ANNOTATED": 3.5, "ANNOTATED_INT": 7}, found
    dropped = {n for _m, n, _r in scan.excluded}
    assert dropped == {"LEDGER_PATH", "NOT_A_NUMBER"}, dropped


def test_the_modules_the_hand_count_missed_are_all_present() -> None:
    """The six the hand inventory never reached, named individually so a regression says which.

    ⚠ These are anchors, not a total. Pinning the total would make every constant added to the tree a
    failure of this test rather than a fact about the tree, which is the opposite of a denominator.
    """
    scan = sweep()
    modules = {r.module for r in scan.rows}
    for stem in ("fa_maturation", "ecm_library", "hand_kmc",
                 "polarization_activegel", "polymerization_warp", "piezo"):
        assert any(m.endswith(f"laws/{stem}.py") for m in modules), stem


def test_no_value_is_produced_by_importing_the_module_under_scan() -> None:
    """Values come from literals only. Importing to read a constant runs the module's import side.

    In this tree that would mean `wp.init()` and a device query from a counting script, so a value the
    parser cannot read is recorded as unreadable rather than evaluated. This pins the refusal.
    """
    scan = sweep()
    for row in scan.rows:
        assert row.value is None or isinstance(row.value, (int, float)), (row.module, row.name)
    assert any(r.value is None for r in scan.rows), (
        "not one declaration was unreadable, which means the parser is evaluating something")


def test_the_roots_are_directories_and_the_sweep_reports_an_empty_one() -> None:
    """A root that yields nothing must be REPORTED, never silently contribute zero."""
    assert ROOTS, "no roots"
    scan = sweep(roots=ROOTS + ("aleph/does_not_exist",))
    assert any("does_not_exist" in note for note in scan.empty_packages), scan.empty_packages


def test_the_kind_hint_reads_the_name_and_claims_nothing_more() -> None:
    """`kind` splits machine knobs from cell physics BY NAME, and both classes must be non-empty.

    ⚠ A hint, not a verdict. If it ever collapses to one class it has stopped distinguishing anything
    and a reader would take `DEFAULT_MIN_TAU_WINDOWS` and `F_HEAD_PN` for the same kind of thing.
    """
    scan = sweep()
    machine = sum(1 for r in scan.rows if r.kind == "machine-hint")
    assert 0 < machine < len(scan.rows), (machine, len(scan.rows))


def test_the_oracle_tree_is_not_swept() -> None:
    """`virtual_cell/` is the RULER and counting it inside N corrupts the question N answers.

    ⚠ PI ruling 2026-08-24, settled by measurement rather than argument: nothing outside
    `virtual_cell/` imports `optics` or `synthetic_microscopy`, and nothing in `world/` / `observe/` /
    `laws/` / `scripts/` imports `virtual_cell` at all. The stage-0 observables are computed from
    arrays geometrically. The optical model enters only for an image-to-image comparison, which is the
    sandbox's question and the external-network data path.
    """
    from aleph.scripts.param_ledger import _NOT_SWEPT
    assert "aleph/virtual_cell" not in ROOTS, ROOTS
    swept_out = {root for root, _reason in _NOT_SWEPT}
    assert "aleph/virtual_cell" in swept_out
    assert all(reason.strip() for _root, reason in _NOT_SWEPT), "a root excluded with no reason"
    scan = sweep()
    assert not any("virtual_cell" in r.module for r in scan.rows), "the oracle tree leaked in"


def test_the_kind_hint_is_documented_as_unfit_for_filtering() -> None:
    """The hint's measured false-positive rate must travel with it, not live in a commit message.

    ⚠ 13 of 21 flagged rows on the physics tree are physical, including `E_SUB_DEFAULT_PA` — a
    stage-0 observable's independent variable, flagged only for containing the word DEFAULT. A reader
    who filters on `kind` drops it. The warning is pinned here so removing it fails.
    """
    import inspect

    from aleph.scripts.param_ledger import as_yaml
    import aleph.scripts.param_ledger as mod

    # ⚠ The measured rate lives in a `#:` comment beside `_MACHINE_HINTS`, so it is read from the
    # SOURCE. An earlier draft of this line was `assert ... if False else True`, which is a check that
    # cannot fail -- the exact defect class this repository keeps finding, written by the person
    # writing the test for it.
    source = inspect.getsource(mod)
    assert "13 of 21" in source or "WRONG 13 times in 21" in source, (
        "the hint's measured false-positive rate left the source; a bare 'this is only a hint' "
        "does not tell a reader it is wrong more often than right")
    text = as_yaml(sweep())
    assert "MAY NOT be used as a filter" in text
    assert "E_SUB_DEFAULT_PA" in text, "the worked example of the hint being wrong went missing"


def test_every_row_carries_an_unassigned_partition() -> None:
    """Decision 3 gets a slot per row, and it starts empty. It sets the real K x M x N total."""
    from aleph.scripts.param_ledger import as_yaml
    text = as_yaml(sweep())
    assert "partition: UNASSIGNED" in text
    for label in ("SHARED", "PER_CELL_LINE", "PER_STATE", "INSTRUMENT"):
        assert f"partition: {label}" not in text, f"the ledger assigned {label} to something"


def test_every_row_lands_ungraded() -> None:
    """The ledger asserts no provenance. Grading is the PI's and a guessed grade is false coverage."""
    from aleph.scripts.param_ledger import as_yaml
    text = as_yaml(sweep())
    assert "grade: UNGRADED" in text
    for grade in ("SOURCED", "DERIVED", "CONVENIENCE", "PI_GAP", "UNRATIFIED_PROXY"):
        assert f"grade: {grade}" not in text, f"the ledger assigned {grade} to something"


def test_the_generated_yaml_parses() -> None:
    """It is hand-formatted for review, so its syntax is checked rather than assumed."""
    import yaml
    from aleph.scripts.param_ledger import as_yaml
    doc = yaml.safe_load(as_yaml(sweep()))
    assert doc["n_declarations"] == len(doc["rows"])
    assert doc["n_excluded"] == len(doc["excluded"])
    assert ast.literal_eval(repr(doc["rows"][0]["line"])) > 0
